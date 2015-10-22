#!/usr/bin/env perl
# vcf2pseudoalignmnet
# Purpose:  Given a set of *.vcf files, examines all called SNPs and generates a 'pseudoalignment'
#  of all the high quality SNPs/base calls (to be used for further phylogenetic analysis).

use warnings;
use strict;
use FindBin;
use lib $FindBin::Bin.'/lib';
use lib $FindBin::Bin;
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;
use Streaming;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::SimpleAlign;
# VCF module from vcftools: http://vcftools.sourceforge.net/index.html
use Vcf;
use File::Temp qw /tempdir/;
use InvalidPositions;
use List::MoreUtils qw/all any firstidx/;
use File::Path qw /rmtree /;



my $verbose;

sub usage
{
	"Usage: $0 --vcf-dir [vcf dir] --mpileup-dir [mpileup dir] --output-base [base of output alignment file]\n".
	"Parameters:\n".
	"\t--vcf-dir: The directory containing the vcf files.\n".
        "\t--vcfsplit: Multiple list of key/value pair  'name=path/to/vcf.gz'".          
	"\t--mpileup-dir: Directory containing the vcf files produced by 'samtools mpileup \$file | bcftools view -cg' (*.vcf).\n".
        "\t--mpileup: Multiple list of key/value pair  'name=path/to/vcf.gz'".
	"\t-o|--output-base:  The output base name for the alignment file(s)\n".
	"Options:\n".
	"\t-r|--reference:  The name of the reference to use in the alignment (default: reference)\n".
	"\t-f|--format:  The format to output the alignment to, one of the Bio::AlignIO supported formats (default: fasta)\n".
	"\t-c|--coverage-cutoff:  The cutoff for coverage to include a reference base (default: 1)\n".
	"\t--invalid-pos: A TSV file that contains a list of range(s) (one per line) of CHROM\\tSTART_POS\\tEND_POS\\n".
	"\t--verbose:  More information printed\n".
	"\t-h|--help:  Help\n";
}

sub create_mpileup_table
{
	my ($vcf_files, $vcf_dir, $mpileup_dir) = @_;
	my %mpileup_table;

	for my $vcf_name (keys %$vcf_files)
	{
		my $vcf_file = $vcf_files->{$vcf_name};
		my $mpileup_file = "$mpileup_dir/".basename($vcf_file);
		$mpileup_table{$vcf_name} = $mpileup_file;
		if (not (-e $mpileup_file))
		{
			print STDERR "Could not find mpileup file \"$mpileup_file\" corresponding to vcf-file \"$vcf_file\"\n";
			return undef;
		}
	}

	return \%mpileup_table;
}


############
### MAIN ###
############

# maps format name to format file extension
my %valid_formats = ('fasta' => 'fasta', 'phylip' => 'phy', 'clustalw' => 'cl');


my ($vcf_files,$mpileup_files,$coverage_cutoff,$bcftools,$requested_cpus,$output_base,$formats,
    $refs_info,$invalid_pos,$invalid_total,$reference, $min_mean_mapping, $ao
) = prepare_inputs();

#create temp working directory for all combines vcf files
#in future make them stay around...
my $tmp_dir = tempdir (CLEANUP => 1);
    
#combine the mpileup and freebayes vcf files together
#in the future, might be taken out to it's own script
my $files = combine_vcfs($vcf_files,$mpileup_files, $coverage_cutoff,$bcftools,$tmp_dir,$requested_cpus,$min_mean_mapping,$ao);


my $valid_positions = $output_base . "-positions.tsv";


#create pseudo-positions.tsv file
my $stats = filter_positions($files,$refs_info,$invalid_pos,$valid_positions,$bcftools,$requested_cpus);

#create pseudo-stats.csv file
print_stats($stats,$invalid_total,$output_base . '-stats.csv');


#create alignment files
for my $format (@{$formats})
{
    my $output_file = $output_base . ".".$valid_formats{$format};
    print STDERR "Alignment written to $output_file\n";
    my $cmd = "$FindBin::Bin/positions2pseudoalignment.pl -i $valid_positions -f $format --reference-name $reference -o $output_file";
    print "$cmd\n";
    system($cmd) == 0 or die "Could not run $cmd";
}



exit;


sub combine_vcfs{
    my ($vcf_files,$mpileup_files, $coverage_cutoff,$bcftools,$tmp_dir,$cpus,$min_mean_mapping,$ao) = @_;

    my %files;

    my $pm;
    my $num_cpus=`cat /proc/cpuinfo | grep processor | wc -l`;
    chomp $num_cpus;
    #ensure that you user cannot request more threads then CPU on the machine
    if ( $cpus > $num_cpus) {
        $cpus = $num_cpus;
    }

    $pm=Parallel::ForkManager->new($cpus);
    
    $pm -> run_on_finish ( # called BEFORE the first call to start()
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $child_data) = @_;
            # retrieve data structure from child
            if (defined($child_data)) {  # children are forced to send anything
                my ($name) = keys %$child_data;
                $files{$name} = $child_data->{$name};
            } else {
                die "One or more vcf file did not produce any data!\n";
            }
        }
    );

    #intermediate files will be in vcf or bcf format. Adding the ability to hardcode the switch because would like to keep using bcf because of space and speed but having soooo much trouble with issues that
    #we need to switch to using vcf. Hence, will add ability to toggle between the two formats with a single commmenting one line
    #issue has been reported on github bcftools as issue # 317
#   my ($ext,$out_type) = ('.bcf','-O b');
    my ($ext,$out_type) = ('.vcf.gz','-O z');

    
    foreach my $sample( keys %$mpileup_files) {
        
        my $pid = $pm->start and next;
        my $f_file = $vcf_files->{$sample};
        my $m_file = $mpileup_files->{$sample};

        my $cmd;
        
        my $file_name = "$tmp_dir/$sample" . "_combined.bcf.gz";

        my ($dir) = "$tmp_dir/$sample" . '_answer';

        if ( not -d $dir) {
            mkdir $dir;
        }

        

        #before running anything else
        #need to run filtered-coverage on original mpileup bcf file
        #going to use the default one provided from bcftools filter instead of a custom one.
        #if we are using vcf files, need to convert first to vcf before applying the filter. The reason is bug with bcftools where SOME files will put the wrong FLAG in....
        if ( $ext eq '.vcf.gz') {
            $cmd = "$bcftools view $m_file -O v | $bcftools  filter -s 'coverage' -i 'DP>=$coverage_cutoff'  $out_type > $dir/coverage_mpileup$ext";
        }
        else {
            $cmd = "$bcftools  filter -s 'coverage' -i 'DP>=$coverage_cutoff' $m_file $out_type > $dir/coverage_mpileup$ext";            
        }

        system($cmd) == 0 or die "Could not run $cmd";


        $cmd = "$bcftools index  $dir/coverage_mpileup$ext";
        system($cmd) == 0 or die "Could not run $cmd";

        
        #confirm SNPS in freebayes by comparing them to mpileup REF and ALT columns
        #mark all SNPS found in mpileup but NOT in freebayes as filtered-mpileup with 'some' option. 'some' options allows
        #only records where some subset of ALT alleles match are compatible
        #                  #chrom      #pos     #ref  #alt 
        #so if mpileup had NC_007530.2|668709 . T     G,A
        #it will map to a freebayes with
        #                  NC_007530.2|668709 . T     G
        $cmd = "$bcftools  isec $f_file $dir/coverage_mpileup$ext -p $dir -c some $out_type";
        system($cmd) == 0 or die "Could not run $cmd";

	#result from bcftools isec is a directory that contains multiple *$ext files
	#ignoring: 0000$ext is for unique position just for vcf-split aka freebayes (we do not use at all since should be empty because mpileup report ALL positions)
	#use: 0001$ext is for unique position only found in mpileup 
	#use: 0002$ext is positions that both freebayes and mpileup have a consensus on the base pair call (either a SNP or same as reference)
	      #VCF line that is kept is the one from freebayes and NOT mpileup
	#use: 0003$ext same as 0002$ext but where mpileup VCF line is kept and not freebayes. Need so we can confirm isec SNPS from freebayes (0002$ext)


        #filter out SNPs that were only found in mpileup but NOT in freebayes
        $cmd = "$bcftools  filter  -m + -s 'mpileup' -i ' TYPE!=\"snp\" ' $dir/0001$ext $out_type  > $dir/1-0001$ext";
        
        system($cmd) == 0 or die "Could not run $cmd";


	


	#need to get rid of the stupid format information since they cannot be merge later downstream
	$cmd ="$bcftools  view -h $dir/1-0001$ext";
        my $result = `$cmd`;
	if ($result =~ /##FORMAT=\<ID=GL/){
	    $cmd = "$bcftools  annotate -x FORMAT -x FORMAT/GT -x FORMAT/GL  $dir/1-0001$ext $out_type > $dir/filtered_mpileup$ext";
	}
        elsif ( $result eq '') {
            die "Failed to retrieve header of '1-0001$ext' for strain '$sample'\n";
        }
	else{
	    $cmd = "$bcftools  annotate -x FORMAT -x FORMAT/GT $dir/1-0001$ext $out_type > $dir/filtered_mpileup$ext";
	}
        system($cmd) == 0 or die "Could not run $cmd";

	######################################################################################################

       
        #Doing two level of filtering here
        #first is by alternative allele ratio. Default is 75% of the reads need to show a SNP exist
        #Second is MQM (min mean mapping) needs to be above a threshold (default is 30)
        #if either of them fail, it will be hard clip out.
        #NB that not sure how it handles when have multiple different alternative alleles
        #also hard clipping ones that fail filtering. Do not want to have them appear in the pseudo-positions since they never passed
        $cmd = "$bcftools  filter  -m + -e  'MQM<$min_mean_mapping || AO/DP<$ao'  $dir/0002$ext $out_type   > $dir/1-0002$ext && bcftools index $dir/1-0002$ext";
        system($cmd) == 0 or die "Could not run $cmd";



        my $mpileup_checked_bcf = check_reference($bcftools,"$dir/1-0002$ext","$dir/0003$ext",$dir,"$dir/filtered_freebayes$ext",$out_type);

        if ($mpileup_checked_bcf ) {
            die "Could not corretly format intersection mpileup file\n";
        }


	$cmd = "$bcftools  annotate -x FORMAT -x FORMAT/GL -x FORMAT/GQ  $dir/filtered_freebayes$ext $out_type > $dir/filtered_freebayes2$ext";
        system($cmd) == 0 or die "Could not run $cmd";

        
        $cmd = "$bcftools index  $dir/filtered_freebayes2$ext";
        system($cmd) == 0 or die "Could not run $cmd";


        $cmd = "$bcftools index  $dir/filtered_mpileup$ext";
        system($cmd) == 0 or die "Could not run $cmd";

        #need to have mpileup header otherwise bcftools index has issues
        #reason being is that freebayes does not report ##contig which is needed for bcftools index/query. Issue arise with some edge cases
        $cmd = "$bcftools merge --print-header  $dir/filtered_mpileup$ext $dir/filtered_freebayes2$ext > $dir/header";
        system($cmd) == 0 or die "Could not run $cmd";

        $cmd = "$bcftools  merge -O b --use-header $dir/header $dir/filtered_freebayes2$ext $dir/filtered_mpileup$ext > $file_name";
        system($cmd) == 0 or die "Could not run $cmd";

        $cmd = "$bcftools index -f $file_name";
        system($cmd) == 0 or die "Could not run $cmd";
        
        rmtree $dir;
        $pm->finish(0,{"$sample" =>$file_name});        
    }
    
    $pm->wait_all_children;

    
    return \%files;
}


sub check_reference {
    my ($bcftools,$freebayes,$mpileup,$basedir,$output,$out_type) = @_;
    my $cmd;
    
        #check to see if we have any records to run again
        my $stats = `$bcftools  stats  $freebayes`;
        if ( $stats =~ /number of records:\s+(\d+)/) {
            if ($1) {
                #check to ensure that there is no reference in the header that does NOT at least one have record. It there is no record for a reference
                #bcftools will simply freeze up
                #get header line and parse out the ##contig=

                my $out = `$bcftools view -h $mpileup | grep "##contig" | sed -e 's/^##contig=<ID=//' -e 's/,length=[0-9][0-9]*>//'`;
                die "Error: no ##contig entries in '$out'" if ($out =~ //); #assume that we have at least one reference
                
                my %refs;
                
                foreach ( split /\n/,$out) {
                    $refs{$_}++;
                }
                
                $out = `$bcftools index -s  $mpileup`;
                die "Error: no entry found using bcftools index  in '$out'" if ($out =~ //); #assume that we have at least one reference coming back from bcftools index
                foreach my $line( split/\n/,$out) {
                    my @data=split/\t/,$line;
                    if (exists $refs{$data[0]}) {
                        delete $refs{$data[0]};
                    }
                }

                #if any reference still in %refs, means there is NO position in the vcf and need to be removed!
                if ( keys %refs) {
                    my $old_header = `$bcftools view -h $mpileup`;
                    open my $out, ">$basedir/newheader";

                    foreach my $line( split/\n/,$old_header) {
                        my $filter_in=1;
                        foreach my $ref( keys %refs) {
                            if ( (index $line,$ref) !=-1) {
                                $filter_in=0;
                                last;
                            }
                        }
                        if ( $filter_in) {
                            print $out "$line\n";
                        }

                    }
                    close $out;
                    #put new header on the file
                    $cmd = "$bcftools view -H  $mpileup > $basedir/beer";
                    system($cmd) == 0 or die "Could not run $cmd";
                    $cmd = "cat $basedir/newheader $basedir/beer > $basedir/beer2";
                    system($cmd) == 0 or die "Could not run $cmd";
                    $cmd = "$bcftools view $basedir/beer2 $out_type > $mpileup";
                    system($cmd) == 0 or die "Could not run $cmd";
                    $cmd = "$bcftools index -f $mpileup";
                    system($cmd) == 0 or die "Could not run $cmd";
                }
                
                $cmd = "$bcftools  annotate  $freebayes -a $mpileup $out_type -c FILTER     > $output";
                system($cmd) == 0 or die "Could not run $cmd";
            }
            else {
                $cmd = "ln -s $freebayes  $output";
                system($cmd) == 0 or die "Could not run $cmd";
            }
    }
    else {
        die "Could not run bcftools stats on file '$freebayes'\n";
    }

    
        

        
    return 0;
}


sub filter_positions {
    my ($files,$refs,$invalid_pos,$valid_positions,$bcftools,$cpus) = @_;

    my %results;
    my $pm;
    my $num_cpus=`cat /proc/cpuinfo | grep processor | wc -l`;
    chomp $num_cpus;
    #ensure that you user cannot request more threads then CPU on the machine
    if ( $cpus > $num_cpus) {
        $cpus = $num_cpus;
    }

    my %vcfcore;
    

    $pm=Parallel::ForkManager->new($cpus);
    
    $pm->run_on_finish ( # called BEFORE the first call to start()
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $child_data) = @_;
            # retrieve data structure from child
            if (defined($child_data)) {  # children are forced to send anything
                
                $results{$child_data->{'job_file'}} = $child_data->{'job_file'}; # get file_name

                #get stats information
                foreach my $chrom ( keys %{$child_data->{'stats'}} ) {
                    $vcfcore{$chrom}{'invalid'} += $child_data->{'stats'}{$chrom}{'invalid'};
                    $vcfcore{$chrom}{'core'}    += $child_data->{'stats'}{$chrom}{'core'};
                }
            } else {
                die "One or more vcf file did not produce any data!\n";
            }
        }
    );
    my $parallel=1;
    
    my $bit_size = 100000;
    my $job_id=0;
    my $tmp_dir = tempdir (CLEANUP => 1);

    #print header file
    open my $out, '>', "$tmp_dir/$job_id";
    print $out "#Chromosome\tPosition\tStatus\tReference\t";
    my @samples_list = sort {$a cmp $b } keys %$files;
    print $out join("\t",@samples_list);
    print $out "\n";
    close $out;
    
    
    foreach my $chrom( sort { $a cmp $b } keys %$refs) {
        
        my ($start,$stop)=(0,0);
        my ($range);
        my $length = $refs->{$chrom}{'length'};
        $vcfcore{$chrom}= {'invalid' => 0,
                           'core' => 0,
                           'total'=>$length
                       };
        
        if ( !$length) {
            die "Could not determine length of '$chrom' from provided reference fasta file\n";
        }
        
        
        while ($stop < $length ) {
            #inclusive range, we do not care if we go over b/c vcf query will only return what it can.
            $start = $stop +1;
            $stop = $stop+$bit_size;
            $range="$chrom:" . join('-',($start,$stop));

            $range = quotemeta $range;

            $job_id++;

            my $pid = $pm->start and next if $parallel;
            
            my %stats = ($chrom => { 'invalid' => 0,
                                    'core' => 0,
                                });
            
            my $f_name = "$tmp_dir/$job_id";
            open my $out, '>',$f_name;
            
            my $streamers = Streaming::create_streamers($files,$range,$job_id,$bcftools);
            my $cur_pos = $start;
            
            my @data = $streamers->($chrom,$cur_pos);


            #search for all entries to have "EOF" as their data.
            while (any { $_->{'status'} ne 'EOF' } @data) {

                #get current bp for current positions
                my $ref_bp = $refs->{$chrom}{'bps'}[$cur_pos];

                #if we do not have any SNP, we do not print at all!
                if ( any {  exists $_->{'alt'} && $_->{'alt'} ne '.' } @data ) {
                    my @line = ($chrom,$cur_pos);

                    #all columns are passed all cut off parameters
                    if ( all { $_->{'status'} eq 'PASS' }  @data) {

                        if ($invalid_pos && exists $invalid_pos->{"${chrom}_${cur_pos}"} ) {
                            push @line, 'filtered-invalid';
                            $stats{$chrom}{'invalid'}++;
                        }
                        else {
                            push @line, 'valid';
                            $stats{$chrom}{'core'}++;
                        }
                        
                        push @line,$ref_bp;
                        foreach ( @data) {
                            if ( $_->{'alt'} eq '.') {
                                push @line,$ref_bp;
                            }
                            else {
                                push @line,$_->{'alt'};
                            }
                        }
                        print $out join("\t",@line) . "\n";
                        
                    }
                    #now we are working with position where there is at least one SNP and at least one filtered-*
                    elsif (any {  ( $_->{'status'} eq 'PASS' || $_->{'status'} eq 'filtered-invalid' )  && $_->{'alt'} ne '.'} @data ){
                        
                        my $is_core=1; #if all positions have coverage, part of the core regardless if we do not agree it's a high quality SNP or invalid position                
                        if ($invalid_pos && exists $invalid_pos->{"${chrom}_${cur_pos}"} ) {
                            push @line, 'filtered-invalid';
                            $stats{$chrom}{'invalid'}++;
                            $is_core=0;
                        }
                        else {
                            #always take filtered-invalid before anything else if it there.
                            my $index = firstidx {$_->{'status'} eq 'filtered-invalid'  } @data;

                            #if we did not find a filtered-invalid, find another none valid status
                            if ( $index == -1) {
                                $index = firstidx {$_->{'status'} ne 'PASS' } @data;
                            }
                            
                            if ( $data[$index]->{'status'} eq '' || $data[$index]->{'status'} eq 'EOF') {
                                push @line,'filtered-coverage';
                            }
                            else {
                                my $stats =  $data[$index]->{'status'};
                                $stats = 'filtered-' . $stats;
                                push @line,$stats;
                            }
                        }

                        
                        push @line,$ref_bp;

                        
                        foreach my $col ( @data) {
                            my $status = $col->{'status'};
                            
                            #if '', implies there is no reads covering that $cur_pos
                            if ( $status eq '' || $status eq 'EOF' || $status eq 'coverage') {
                                push @line,'-';
                                $is_core=0;
                            }
                            else {
                                #if we have filtered-mpileup, we have inconsistent calles between variant callers
                                if ($status eq 'mpileup' ) {
                                    push @line,'N';
                                }
                                #have a position that passes the cut-off parameter. Either show the SNP or the reference
                                elsif ( $status eq 'PASS' || $status eq 'filtered-invalid') {
                                    $is_core=0 if $status eq 'filtered-invalid';
                                    
                                    if ( $col->{'alt'} eq '.') {
                                        push @line,$ref_bp;
                                    }
                                    else {
                                        push @line,$col->{'alt'};
                                    }
                                }
                                else {
                                    push @line,'-';
                                }
                                
                            }
                        }
                        
                        $stats{$chrom}{'core'}++ if $is_core; #if all positions had coverage, add to the core
                        
                        print $out join("\t",@line) . "\n";
                        
                        
                    }#end else
                    else {
                        print "WARNING: Edge case found at '$chrom' on position '$cur_pos'\n" if $verbose;
                    }
                    
                }#end if
                #see if everyone has a PASS status and there is no SNP
                elsif ( all { $_->{'status'} eq 'PASS' }  @data ) {

                    #check to see if the positions in the invalid position, if not, count it toward the total
                    if (not ($invalid_pos && exists $invalid_pos->{"${chrom}_${cur_pos}"} )) {
                        $stats{$chrom}{'core'}++;                                            
                    }

                    
                }
                else {
                    print "WARNING: Edge case found at '$chrom' on position '$cur_pos'\n" if $verbose;
                }
                
                
                $cur_pos++;
                @data = $streamers->($chrom,$cur_pos);
            }
            
            close $out;
            
            $pm->finish(0,{'job_file' =>$job_id, 'stats' => \%stats}) if $parallel;
            $results{$job_id}=$job_id if not $parallel;
        }
        
    
    }
    $pm->wait_all_children if $parallel;   

    #combine all results files in order
    my @cmd = "cat $tmp_dir/0";
    foreach ( sort {$a <=> $b } keys %results) {
         push @cmd, "$tmp_dir/" . $results{$_};
    }
    my $cmd = join(' ' , @cmd) . " > $valid_positions";
    
    system($cmd) == 0 or die "Could not run $cmd";
    unlink "0";
    
    map { unlink $_ } values %results;


    
    return \%vcfcore;
}

sub refs_info {
    my ($file) = @_;


    my %refs;
    my $in = Bio::SeqIO->new(-format=>'fasta',-file=>$file);
    while ( my $seq = $in->next_seq()) {
        $refs{$seq->display_id()}{'length'} = $seq->length;
        my @bps = ('X'); #seqs start at index 1 and not zero. So putting a value.
        
        map { push @bps, uc $_ } split //, $seq->seq;
        $refs{$seq->display_id()}{'bps'} = \@bps;
    }


    return \%refs;
    
}


sub print_stats {
    my ($stats,$invalid_total,$out) = @_;

    my %vcfcore = %{$stats};
    
    open my $out_fh,'>',$out;

    my @header = ("#Reference","total length","total invalid pos","total valid pos","total core","Percentage in core");
    
    print $out_fh join("\t",@header) . "\n";
    #negative final_invalid indicates there was no invalid positions file given
    my ($final_core,$final_total,$final_invalid)= (0,0,0);

    if (not  $invalid_total) {
        $final_invalid='N/A';
    }
    
    foreach my $chrom( sort {$a cmp $b } keys %vcfcore) {
        my ($core,$total) = ($vcfcore{$chrom}{'core'},$vcfcore{$chrom}{'total'});
        my $invalid = 'N/A';
        my $total_no_invalid;
        
        my $perc;

        
        #if there was invalid found
        if (exists $invalid_total->{$chrom}) {
            $invalid = $invalid_total->{$chrom}; #grab the total number of invalid positions marked by the user for the current chrom
            $final_invalid +=$invalid;
            
            if ($total == $invalid ) {
                $perc ="NaN";
                $total_no_invalid=0;
            }
            else {
                $perc = sprintf("%.2f",($core/( $total-$invalid)*100));
                $total_no_invalid = $total-$invalid;
            }
            
        }
        elsif ( $vcfcore{$chrom}{'invalid'} && not exists $invalid_total->{$chrom}) {
            die "Found invalid positions for '$chrom' but was not found in invalid positions file.\n";
        }
        else {
            $perc = sprintf("%.2f",($core/$total*100));
            $total_no_invalid = $total;
        }
        $final_core +=$core;
        $final_total +=$total;

        print $out_fh join ("\t", ($chrom,$total,$invalid,$total_no_invalid,$core,$perc)) . "\n";
        
    }


    #getting the total for all references
    my ($perc,$final_total_no_invalid);
    
    if ($final_invalid eq 'N/A' || $final_invalid == 0 ) {
        $perc = sprintf("%.2f",$final_core/$final_total*100);
        $final_total_no_invalid=$final_total;
    }
    else {
        if ($final_total == $final_invalid ) {
            $perc ="NaN";
            $final_total_no_invalid=0;
         }
        else {
            $perc = sprintf("%.2f",$final_core/($final_total - $final_invalid)*100);
            $final_total_no_invalid=$final_total-$final_invalid;
        }
    }
    print $out_fh join ("\t",'all',$final_total,$final_invalid,$final_total_no_invalid,$final_core,$perc) . "\n";

    return;
}



sub prepare_inputs {

    my ($vcf_dir, $mpileup_dir, $output_base, @formats, $reference, $coverage_cutoff, $min_mean_mapping, $ao);
    my ($help, $requested_cpus, $invalid,$fasta,$bcftools,$inc_list,$exc_list);
    my (%vcf_files, %mpileup_files);



    if (!GetOptions('vcf-dir|d=s' => \$vcf_dir,
                    'vcfsplit=s' => \%vcf_files,
                    'mpileup-dir|b=s' => \$mpileup_dir,
                    'mpileup=s' => \%mpileup_files,
                    'format|f=s' => \@formats,
                    'output-base|o=s' => \$output_base,
                    'reference|r=s' => \$reference,
                    'fasta=s' => \$fasta,
                    'coverage-cutoff|c=i' => \$coverage_cutoff,
                    'invalid-pos=s' => \$invalid,
                    'min-mean-mapping=i'=> \$min_mean_mapping,
                    'ao=s' => \$ao,
                    'help|h' => \$help,
                    'numcpus=i' => \$requested_cpus,
                    'b|bcftools-path=s' => \$bcftools,
                    'verbose|v' => \$verbose,
                    'include=s' => \$inc_list,
                    'exclude=s' => \$exc_list)){
        die "Invalid option\n".usage;
    }
    
    print usage and exit(0) if (defined $help);
    $verbose = 0 if (not defined $verbose);

    
    if ( $vcf_dir and $mpileup_dir){
        die "vcf-dir does not exist\n".usage if (not -e $vcf_dir);
        die "mpileup-dir does not exist\n".usage if (not -e $mpileup_dir);
    }
    elsif ( scalar keys %vcf_files == 0 or scalar keys %mpileup_files ==0){
        die "Was not able to find any vcf files from freebayes and/or mpileup.";
    }
    
    die "output-base undefined\n".usage if (not defined $output_base);
    
    if (not defined $bcftools or not -e $bcftools){
        
        #check to see if bcftools is on the path
        #normally we always want to be passed the path to the tool but this is for Galaxy implementation
        my $alive=`bcftools 2>&1 1>/dev/null`;
        if ( $alive && $alive =~ /Program: bcftools/) {
            $bcftools="bcftools";
        }
        else {
            die "bcftools-path not defined and not found on the path.\n".usage         
        }
        
    }
    
    if(not defined $min_mean_mapping){
    	$min_mean_mapping = 30;
    }

    #need to have ratio/percentage in double format i.e 75% = 0.75
    if(not defined $ao){
       $ao = 0.75;
    }
    elsif ( $ao > 1) {
        print "Assuming that '$ao' is given as percentage for alternative allele. Changing to a decimal\n";
        $ao = $ao/100;
    }
    
    #need check to see if bcftools was complied with htslib and also has the correct plugin installed
    my $usage_state = `$bcftools 2>&1 1>/dev/null`;
    if ( not $usage_state =~ /Version: .* \(using htslib/ ) {
        die "bctools was not complied with htslib.\nPlease re-compile with htslib\nInstruction: http://samtools.github.io/bcftools/\n";
    }
    
    if (not defined $reference){
        print STDERR "reference name not defined, calling it 'reference'\n";
        $reference = 'reference';
    }

    
    
    if (defined $invalid){
        if ( ! -e $invalid){
            die "Was given an invalid position file but could not locate it '$invalid'\n";
        }
    }
    else{
        print STDERR "invalid position file not defined, Will ignore step\n";
    }
    
    
    $requested_cpus = 1 if (not defined $requested_cpus);
    
    if (@formats <= 0){
	print STDERR "warning: format not defined, assuming fasta\n";
	@formats = ("fasta");
    }
    else{
        for my $format (@formats){
            die "unrecognized format '$format', must be one of '".join(' ', keys %valid_formats),"'\n" if (not defined $valid_formats{$format});
        }
    }
    
    if (not defined $coverage_cutoff){
        print STDERR "warning: coverage-cutoff not set, assuming it is 1\n";
        $coverage_cutoff = 1;
    }
    elsif ($coverage_cutoff !~ /^\d+$/){
        die "coverage-cutoff=$coverage_cutoff is invalid\n".usage;
    }
    
    
    my $dh;
    # fill table vcf_files with entries like if only provided the vcf-dir
    #  vcf1 => dir/vcf1.bcf.gz
    #  vcf2 => dir/vcf2.bcf.gz
    if ( $vcf_dir) {
        opendir($dh, $vcf_dir) or die "error opening directory $vcf_dir: $!";
        %vcf_files = map { /^(.*)\.bcf\.gz$/; $1 => "$vcf_dir/$_"} grep { /\.bcf\.gz$/ } readdir($dh);
        closedir($dh);
    }
    
    die "No *.bcf.gz files found in $vcf_dir.  Perhas you need to compress and index with 'tabix' tools\n".
        "Example: bgzip file.vcf; tabix -p vcf file.bcf.gz" if (keys(%vcf_files) <= 0);
    
    my $total_samples = (keys %vcf_files);


    if ( $mpileup_dir){
        # create table of mpileup files corresponding to input freebayes/variant vcf files
        # assumes files are named with the same prefix
        # ex vcf-dir/file1.bcf.gz and mpileup-dir/file1.bcf.gz    
        my $mpileup_table = create_mpileup_table(\%vcf_files, $vcf_dir, $mpileup_dir);
        
        if (not defined $mpileup_table){
            die "Error: vcf-dir contains unmatched files in mpileup-dir";
        }
        else{
            %mpileup_files = %{$mpileup_table};
        }
        
    }
    else{
        if (scalar keys %mpileup_files != scalar keys %vcf_files ){
            die "Error: vcfsplit contains uneven number compare to mpileup-files";
        }
        my $m_name= join ('',sort {$a cmp $b } keys %mpileup_files);
        my $v_name= join ('',sort {$a cmp $b } keys %vcf_files);
        if ( $m_name ne $v_name) {
            die "Error: vcfsplit contains unmatched files to mpileup-files";
        }
            
    }
    if ( not -e $fasta) {
        die "Error: Was not given reference fasta file\n";
    }
    
    my $refs_info = refs_info($fasta);
    
    my ($invalid_pos,$invalid_total);
    
    if ($invalid){
        my $invalid_positions_parser = InvalidPositions->new;
        ($invalid_pos,$invalid_total) = $invalid_positions_parser->read_invalid_positions($invalid);
    }
    if ( $inc_list && $exc_list) {
        die "Cannot have both an include and exclude list. Only one please\n";
    }

    #copying reference to current hashes of vcf files
    my ($vcf_files,$mpileup_files)= (\%vcf_files,\%mpileup_files);

    #if we have a including list, then use it
    if ( $inc_list) {
        ($vcf_files,$mpileup_files) = strain_selection($inc_list,1,\%vcf_files,\%mpileup_files);
    }

    #if we have a excluding list, then use it
    if ( $exc_list) {
        ($vcf_files,$mpileup_files) = strain_selection($exc_list,0,\%vcf_files,\%mpileup_files);
    }
    

    return ($vcf_files,$mpileup_files,$coverage_cutoff,$bcftools,$requested_cpus,$output_base,\@formats,
            $refs_info,$invalid_pos,$invalid_total,$reference, $min_mean_mapping, $ao);
}

sub strain_selection {
    my ($list,$keep,$vcf_files,$mpileup_files)=@_;

    my $vcfs = $vcf_files;
    my $mpileups = $mpileup_files;
    my (%tokeep_vcf,%tokeep_mpileup);

    open my $fh, '<',$list || die "Could not open file '$list'\n";
    while (my $name = <$fh>) {
        chomp $name;

        #check to see if the name is in the list of keys
        if ( exists $vcfs->{$name} && exists $mpileups->{$name}) {

            #if we are doing a inclusion list, add it to the new hashes
            #otherwise just delete from the hash coming in
            if ( $keep) {
                $tokeep_vcf{$name} = $vcf_files->{$name};
                $tokeep_mpileup{$name} = $mpileup_files->{$name};
            }
            else {
                delete $vcfs->{$name};
                delete $mpileups->{$name};
            }
        }
        else {
            print STDERR "WARNING: Could not find strain '$name' in list of strains. Ignorning\n";
        }
    }
    close $fh;

    #if we are doing an eclusion list, then move return hashes
    if (not $keep) {
        return ($vcfs,$mpileups);
    }
    else {
        return (\%tokeep_vcf,\%tokeep_mpileup);
    }
    

}

