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

sub variant_info_to_hash
{
	my ($info_string) = @_;

	my @info = split(/;/,$info_string);
	my %info_hash = ();
	foreach my $info_entry (@info)
	{
		if (defined $info_entry and $info_entry ne '')
		{
			my @parts = split(/=/, $info_entry);
			if (defined $parts[0] and $parts[0] ne '')
			{
				if (@parts > 1)
				{
					$info_hash{$parts[0]} = $parts[1];
				}
				else
				{
					$info_hash{$parts[0]} = undef;
				}
			}
		}
	}

	return \%info_hash;
}



############
### MAIN ###
############

# maps format name to format file extension
my %valid_formats = ('fasta' => 'fasta', 'phylip' => 'phy', 'clustalw' => 'cl');


my ($vcf_files,$mpileup_files,$coverage_cutoff,$bcftools,$requested_cpus,$output_base,$formats,
    $refs_info,$invalid_pos,$reference
) = prepare_inputs();

#create temp working directory for all combines vcf files
#in future make them stay around...
my $tmp_dir = tempdir (CLEANUP => 0);
    
#combine the mpileup and freebayes vcf files together
#in the future, might be taken out to it's own script
my $files = combine_vcfs($vcf_files,$mpileup_files, $coverage_cutoff,$bcftools,$tmp_dir,$requested_cpus);


my $valid_positions = $output_base . "-positions.tsv";


#create pseudo-positions.tsv file
my $stats = filter_positions($files,$refs_info,$invalid_pos,$valid_positions,$requested_cpus);

#create pseudo-stats.csv file
print_stats($stats,$output_base . '-stats.csv');


#create alignment files
for my $format (@{$formats})
{
    my $output_file = $output_base . ".".$valid_formats{$format};
    print STDERR "Alignment written to $output_file\n";
    my $cmd = "$FindBin::Bin/positions2pseudoalignment.pl -i $valid_positions -f $format --reference-name $reference -o $output_file";
    print "$cmd\n";
    my $status = system($cmd);
}



exit;


sub combine_vcfs{
    my ($vcf_files,$mpileup_files, $coverage_cutoff,$bcftools,$tmp_dir,$cpus) = @_;

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

    foreach my $sample( keys %$mpileup_files) {
        
        my $pid = $pm->start and next;
        my $f_file = $vcf_files->{$sample};
        my $m_file = $mpileup_files->{$sample};

        my ($cmd,$status);
        
        my $file_name = "$tmp_dir/$sample" . "_combined.vcf.gz";

        my ($dir) = "$tmp_dir/$sample" . '_answer';

        #confirm SNPS in freebayes by comparing them to mpileup REF and ALT columns
        #mark all SNPS found in mpileup but NOT in freebayes as filtered-mpileup with 'some' option. 'some' options allows
        #only records where some subset of ALT alleles match are compatible
        #                  #chrom      #pos     #ref  #alt 
        #so if mpileup had NC_007530.2|668709 . T     G,A
        #it will map to a freebayes with
        #                  NC_007530.2|668709 . T     G
        $cmd = "$bcftools  isec $f_file $m_file -p $dir -c some -O b";
        $status = system($cmd);


        #filter with C complied nml specific filtering, also removing all information in FORMAT column, otherwise we cannot merge farther down
        $cmd = "$bcftools  annotate -x FORMAT -p filter_mpileup:dp=$coverage_cutoff $dir/0001.bcf -O b  > $dir/filtered_mpileup.bcf";
        $status = system($cmd);

        
        #filter by coverage and ratio of 75% with alternative allele
        #also filter by MQM flag = minumum mean mapping quality with > 30
        #NB that not sure how it handles when have multiple different alternative alleles
        #also hard clipping ones that fail filtering. Do not want to have them appear in the pseudo-positions since they never passed
        $cmd = "$bcftools  annotate -x FORMAT -p filter_freebayes:dp=$coverage_cutoff:mqm=30:ao=75  $dir/0002.bcf -O b  > $dir/filtered_freebayes.bcf";
        $status = system($cmd);
        
        

        $cmd = "$bcftools index  $dir/filtered_freebayes.bcf";
        $status = system($cmd);


        $cmd = "$bcftools index  $dir/filtered_mpileup.bcf";
        $status = system($cmd);

        #get the fake merge vcf header
        my $header = $FindBin::Bin . '/fake_vcf_header/header';
        my $bottom_header = $FindBin::Bin . '/fake_vcf_header/bottom_header';

        #need to add header specific reference in the vcf files
        #need a better solution!
        $cmd = "zgrep '##contig' $m_file > $dir/contigs";
        $status = system($cmd);
        $cmd = "cat $header $dir/contigs $bottom_header > $dir/header";
        $status = system($cmd);
        ######
        
        $cmd = "$bcftools  merge -O z  --use-header $dir/header $dir/filtered_freebayes.bcf $dir/filtered_mpileup.bcf > $file_name 2>/dev/null";

        $status = system($cmd);
        
        $cmd = "$bcftools index -t -f $file_name";
        $status = system($cmd);
        
        $files{$sample} = $file_name;
        
        #rmtree $dir;
        $pm->finish(0,{"$sample" =>$file_name});        

    
    
    }
    
    $pm->wait_all_children;

    return \%files;
}


sub filter_positions {
    my ($files,$refs,$invalid_pos,$valid_positions,$cpus) = @_;

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


    #print header file
    open my $out, '>', $job_id;
    print $out "#Chromosome\tPosition\tStatus\tReference\t";
    my @samples_list = sort {$a cmp $b } keys %$files;
    print $out join("\t",@samples_list);
    print $out "\n";
    close $out;
    
    
    foreach my $chrom( keys %$refs) {
        
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
            
            my $f_name = $job_id;
            open my $out, '>',$f_name;
            
            my $streamers = Streaming::create_streamers($files,$range,$job_id);
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
                    else {
                        if ($invalid_pos && exists $invalid_pos->{"${chrom}_${cur_pos}"} ) {
                            push @line, 'filtered-invalid';
                            $stats{$chrom}{'invalid'}++;
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
                                push @line,$data[$index]->{'status'};
                            }
                        }

                        
                        push @line,$ref_bp;
                        
                        foreach my $col ( @data) {
                            my $status = $col->{'status'};
                            
                            #if '', implies there is no reads covering that $cur_pos
                            if ( $status eq '' || $status eq 'EOF') {
                                push @line,'-';
                            }
                            else {
                                #if we have filtered-mpileup, we have inconsistent calles between variant callers
                                if ($status eq 'filtered-mpileup' ) {
                                    push @line,'N';
                                }
                                #have a position that passes the cut-off parameter. Either show the SNP or the reference
                                elsif ( $status eq 'PASS' || $status eq 'filtered-invalid') {
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
                        
                        print $out join("\t",@line) . "\n";
                        
                        
                    }#end else
                    
                    
                }#end if
                #see if everyone has a PASS status and there is no SNP
                elsif ( all { $_->{'status'} eq 'PASS' }  @data ) {
                    $stats{$chrom}{'core'}++;
                }
                
                
                $cur_pos++;
                @data = $streamers->($chrom,$cur_pos);
            }
            
            close $out;
            
            $pm->finish(0,{'job_file' =>$f_name, 'stats' => \%stats}) if $parallel;
            $results{$f_name}=$f_name if not $parallel;
        }
        
    
    }
    $pm->wait_all_children if $parallel;   

    #combine all results files in order
    my @cmd = "cat 0";
    foreach ( sort {$a <=> $b } keys %results) {
         push @cmd, $results{$_};
    }
    my $cmd = join(' ' , @cmd) . " > $valid_positions";
    
    `$cmd`;
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
        my @bps = ('X');
        
        map { push @bps, uc $_ } split //, $seq->seq;
        $refs{$seq->display_id()}{'bps'} = \@bps;
    }


    return \%refs;
    
}


sub print_stats {
    my ($stats,$out) = @_;

    my %vcfcore = %{$stats};
    
    open my $out_fh,'>',$out;
    
    print $out_fh "#Reference,total length,total invalid pos, total core,Percentage in core\n";
    my ($final_core,$final_total,$final_invalid)= (0,0,0);

    foreach my $chrom( keys %vcfcore) {
        my ($core,$total) = ($vcfcore{$chrom}{'core'},$vcfcore{$chrom}{'total'});
        my $invalid = 'N/A';
        my $perc;

        if ( $vcfcore{$chrom}{'invalid'}) {
            $invalid = $vcfcore{$chrom}{'invalid'};
            
            if ($total == $invalid ) {
                $perc =0;
            }
            else {
                $perc = sprintf("%.2f",($core/( $invalid)*100));
            }
            
        }
        else {
            $perc = sprintf("%.2f",($core/$total*100));
        }
        $final_core +=$core;
        $final_total +=$total;
        $final_invalid +=$invalid if $invalid ne 'N/A';
        
            
        print $out_fh join (',', ($chrom,$total,$invalid,$core,$perc)) . "\n";
        
    }


    #getting the total for all references
    my ($perc,$invalid_total)= (0);
    if ($final_invalid) {
        if ($final_total == $final_invalid ) {
            $perc =0;
         }
        else {
            $perc = sprintf("%.2f",$final_core/($final_total - $final_invalid)*100);        
        }
    }
    else {
        $perc = sprintf("%.2f",$final_core/$final_total*100);
    }
    print $out_fh join (',','all',$final_total,$final_invalid,,$final_core,$perc) . "\n";

    return;
}



sub prepare_inputs {

    my ($vcf_dir, $mpileup_dir, $output_base, @formats, $reference, $coverage_cutoff);
    my ($help, $requested_cpus, $invalid,$fasta,$bcftools);
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
                    'help|h' => \$help,
                    'numcpus=i' => \$requested_cpus,
                    'b|bcftools-path=s' => \$bcftools,     
                    'verbose|v' => \$verbose)){
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
    
    
    #need check to see if bcftools was complied with htslib and also has the correct plugin installed
    my $usage_state = `$bcftools 2>&1 1>/dev/null`;
    if ( not $usage_state =~ /Version: .* \(using htslib/ ) {
        die "bctools was not complied with htslib.\nPlease re-compile with htslib\nInstruction: http://samtools.github.io/bcftools/\n";
    }
    my $plugins_state = `$bcftools annotate -l`;
    if ( not ( $plugins_state =~ /-- filter_mpileup --/ && $plugins_state =~ /-- filter_freebayes --/ ) ) {
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
    #  vcf1 => dir/vcf1.vcf.gz
    #  vcf2 => dir/vcf2.vcf.gz
    if ( $vcf_dir) {
        opendir($dh, $vcf_dir) or die "error opening directory $vcf_dir: $!";
        %vcf_files = map { /^(.*)\.vcf\.gz$/; $1 => "$vcf_dir/$_"} grep { /\.vcf\.gz$/ } readdir($dh);
        closedir($dh);
    }
    
    die "No *.vcf.gz files found in $vcf_dir.  Perhas you need to compress and index with 'tabix' tools\n".
        "Example: bgzip file.vcf; tabix -p vcf file.vcf.gz" if (keys(%vcf_files) <= 0);
    
    my $total_samples = (keys %vcf_files);





    if ( $mpileup_dir){
        # create table of mpileup files corresponding to input freebayes/variant vcf files
        # assumes files are named with the same prefix
        # ex vcf-dir/file1.vcf.gz and mpileup-dir/file1.vcf.gz    
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
    
    my $invalid_pos;
    
    if ($invalid){
        my $invalid_positions_parser = InvalidPositions->new;
        $invalid_pos = $invalid_positions_parser->read_invalid_positions($invalid);
    }


    return (\%vcf_files,\%mpileup_files,$coverage_cutoff,$bcftools,$requested_cpus,$output_base,\@formats,
            $refs_info,$invalid_pos,$reference);
}    
