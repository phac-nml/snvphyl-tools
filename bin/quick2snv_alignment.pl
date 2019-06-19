#!/usr/bin/env perl

# quick2snv_alignment
# Purpose:  Given a set of combined freebayes and mpileup files, examines all called SNPs and generates a 'snv_alignment'
#  of all the high quality SNPs/base calls (to be used for further phylogenetic analysis).

use warnings;
use strict;
use FindBin;
use lib $FindBin::Bin.'/../lib';
use lib $FindBin::Bin;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::SimpleAlign;
use File::Temp qw /tempdir/;
use File::Path;
use InvalidPositions;

my %valid_formats = ('fasta' => 'fasta', 'phylip' => 'phy', 'clustalw' => 'cl');
my $verbose;
__PACKAGE__->run unless caller;

sub run
{
	my ($consolidate_vcf,$formats,$output_base,$reference,$invalid_pos,$invalid_total,$requested_cpus,$bcftools,$refs_info) = prepare_inputs(@_);

	my $valid_positions = $output_base . "-positions.tsv";


	#create snvalign-positions.tsv file
	my $stats = filter_positions($consolidate_vcf,$refs_info,$invalid_pos,$valid_positions,$bcftools,$requested_cpus);

	#create snvalign-stats.csv file
	print_stats($stats,$invalid_total,$output_base . '-stats.csv');


	# #create alignment files
	for my $format (@{$formats})
	{
	    my $output_file = $output_base . ".".$valid_formats{$format};
	    print STDERR "Alignment written to $output_file\n";
	    my $cmd = "$FindBin::Bin/positions2snv_alignment.pl -i $valid_positions -f $format --reference-name $reference -o $output_file";
	    print "$cmd\n";
	    system($cmd) == 0 or die "Could not run $cmd";
	}

	return;
}

sub modify_bits {
    my ($bits,$alt,$ref_bp,$bps)=@_;

    #if not given an array reference already, implies we have a fresh one to create
    if ( not $bps) {
        $bps=[()];
    }
    
    my $index=0;
    foreach my $bit( split //,$bits) {
        if ($bit) {
            $bps->[$index]=$alt eq '.' ? $ref_bp : $alt;
        }
        elsif ( not $bps->[$index] ) {
            $bps->[$index]=0;
        }
        $index++;
    }
   
    return $bps;
}
    

sub filter_positions 
{
    my ($consolidate_vcf,$refs,$invalid_pos,$valid_positions,$bcftools,$cpus) = @_;

    my %vcfcore;
    my %results;
    foreach my $chrom( sort { $a cmp $b } keys %$refs) 
    {
        my $length = $refs->{$chrom}{'length'};
        $vcfcore{$chrom}= {'invalid' => 0,
                           'core' => 0,
                           'total'=>$length
                       };
    }    
    
    my $pm=Parallel::ForkManager->new($cpus);

    $pm->run_on_finish ( # called BEFORE the first call to start()
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $child_data) = @_;
            # retrieve data structure from child
            if (defined($child_data)) # children are forced to send anything
            {
                $results{$child_data->{'job_file'}} = $child_data->{'job_file'}; # get file_name
                #get stats information
                foreach my $chrom ( keys %{$child_data->{'stats'}} ) 
                {
                    $vcfcore{$chrom}{'invalid'} += $child_data->{'stats'}{$chrom}{'invalid'};
                    $vcfcore{$chrom}{'core'}    += $child_data->{'stats'}{$chrom}{'core'};
                }
            }
            else
            {
                die "One or more vcf file did not produce any data!\n";
            }
        }
    );
    my $parallel=1;


    my $tmp_dir = tempdir (CLEANUP => 1);

    my $job_id=0;

    #print header file
    open my $out, '>', "$tmp_dir/$job_id";
    print $out "#Chromosome\tPosition\tStatus\tReference\t";
    my @samples_list = sort {$a cmp $b } keys %{$consolidate_vcf};
    print $out join("\t",@samples_list);
    print $out "\n";
    close $out;



    #produce consistent list of filepath for samples
    my $samples_path = join " " , map {$consolidate_vcf->{$_}} @samples_list;

    #determine bin size based on total # of samples so we can ensure to keep memory usage low
    #aiming for less then 1G per thread
    #thoughts that we 10 million bp worth in memory per thread
    my $bin_size = int(10_000_000 / scalar @samples_list);


    foreach my $chrom( sort { $a cmp $b } keys %$refs)
    {
        my ($start,$stop)=(1,0);
        my $length = $refs->{$chrom}{'length'};
        my ($range);

        if ( !$length)
        {
            die "Could not determine length of '$chrom' from provided reference fasta file\n";
        }

        
        while ($stop < $length )
        {

            #inclusive range, we do not care if we go over b/c vcf query will only return what it can.
            $start = $stop;
            $stop = $stop+$bin_size;
            $range="$chrom:" . join('-',($start,$stop));

            $range = quotemeta $range;

            $job_id++;

            my $pid = $pm->start and next if $parallel;
        
            my %stats = ($chrom => {'invalid' => 0,
                                    'core' => 0,
                                });

            my $f_name = "$tmp_dir/$job_id";
            open my $out, '>',$f_name;

            #run bcftool isec on all samples
            my $cmd = "$bcftools isec -r $range -O b $samples_path -n +1 -c snps -f PASS -p $tmp_dir/isec_$job_id";
            system($cmd) == 0 or die "Could not run $cmd";


            #open up sites.txt to determine all SNVs and core position that pass filter values done by consolidate_vcfs.pl  
            open my $in , '<',"$tmp_dir/isec_$job_id/sites.txt";

            #getting first line of file and processing it
            my $line = <$in>;
            chomp $line;

            my ($cur_chrom,$cur_pos,$ref_bp,$alt,$bits) =split(/\t/,$line);


            my $prev_key=$cur_chrom . '_' . $cur_pos;
            my %prev;

            $prev{$prev_key}= {'ref_bp' => $ref_bp,
                               'chrom' =>$cur_chrom,
                               'pos' => $cur_pos,
                               'bp' => modify_bits($bits,$alt,$ref_bp)
                           };
        
            while (my $line= <$in>) {
                chomp $line;
                
                ($cur_chrom,$cur_pos,$ref_bp,$alt,$bits) =split(/\t/,$line);
                
                my $cur_key=$cur_chrom . '_' . $cur_pos;
                
                #determine if we seen all line associated with current chrom and position
                if ($cur_key eq $prev_key) {
                    #add information to current pile
                    $prev{$prev_key}{'bp'}=modify_bits($bits,$alt,$ref_bp,$prev{$prev_key}{'bp'});
                }
                else {
                    #we should have all the information,
                    my @bps=@{$prev{$prev_key}{'bp'}};

                    #check to see if just a core position or a SNV
                    my ($core,$snv)=(1,1);
                    my $prev_bp=$prev{$prev_key}{'ref_bp'};
                    foreach my $bp (@bps) {
                        if ( $bp ne $prev_bp) {
                            #check to see if it is a '0' base pair, if so we are completely disregarding it!
                            if (not $bp ) {
                                $snv=0;
                            }
                            $core=0;
                        }
                    }
                
                    my ($prev_chrom,$prev_pos)=($prev{$prev_key}{'chrom'},$prev{$prev_key}{'pos'});
                    if ($core) {
                        if ($invalid_pos && exists $invalid_pos->{$prev_key} ) {
                            $stats{$prev_chrom}{'invalid'}++;
                        }
                        else{
                            $stats{$prev_chrom}{'core'}++;
                        }
                    }
                    elsif ($snv) {
                        #have a valid snv but has it been marked invalid?
                        if ($invalid_pos && exists $invalid_pos->{$prev_key} ) {
                            $stats{$prev_chrom}{'invalid'}++;
                            print $out join("\t",($prev_chrom,$prev_pos,'filtered-invalid',$prev_bp,@bps)). "\n";
                        }
                        else {
                            $stats{$prev_chrom}{'core'}++;
                            print $out join("\t",($prev_chrom,$prev_pos,'valid',$prev_bp,@bps)). "\n";                    
                        }
                        
                    }
            
                    #empty hash and start gathering next position which may or may not be across multiple lines
                    %prev=();
    
                    $prev{$cur_key}= {'ref_bp' => $ref_bp,
                                      'chrom' =>$cur_chrom,
                                      'pos' => $cur_pos,                              
                                      'bp' => modify_bits($bits,$alt,$ref_bp)                                  
                                  };
                    
                    $prev_key=$cur_key;
                }
            }

            #remove temp isec otherwise we can fill up the hardrive really fast
            rmtree("$tmp_dir/isec_$job_id");
            
            if ( $parallel) {
                $pm->finish(0,{'job_file' =>$job_id, 'stats' => \%stats}) if $parallel;
            }
            else {
                $results{$job_id}=$job_id;

                foreach my $chrom ( keys %stats ) {
                    print "invalid: " . $stats{$chrom}{'invalid'} . "\n";
                    print "core: " . $stats{$chrom}{'core'} . "\n";
                    $vcfcore{$chrom}{'invalid'} += $stats{$chrom}{'invalid'};
                    $vcfcore{$chrom}{'core'}    += $stats{$chrom}{'core'};
                }

            }
            
        }#end while

    } #end foreach
    
    $pm->wait_all_children if $parallel;

    #combine all results files in order
    my @cmd = "cat $tmp_dir/0";
    foreach ( sort {$a <=> $b } keys %results) 
    {
         push @cmd, "$tmp_dir/" . $results{$_};
    }
    my $cmd = join(' ' , @cmd) . " > $valid_positions";

    system($cmd) == 0 or die "Could not run $cmd";
    unlink "0";

    map { unlink $_ } values %results;



    
    return \%vcfcore;
}

sub refs_info 
{
    my ($file) = @_;

    my %refs;
    my $in = Bio::SeqIO->new(-format=>'fasta',-file=>$file);
    while ( my $seq = $in->next_seq()) 
    {
        $refs{$seq->display_id()}{'length'} = $seq->length;
        my @bps = ('X'); #seqs start at index 1 and not zero. So putting a value.

        map { push @bps, uc $_ } split //, $seq->seq;
        $refs{$seq->display_id()}{'bps'} = \@bps;
    }

    return \%refs;
}


sub print_stats 
{
    my ($stats,$invalid_total,$out) = @_;

    my %vcfcore = %{$stats};

    open my $out_fh,'>',$out;

    my @header = ("#Reference name","Total length","Total invalid and excluded positions","Total valid and included positions","Total valid and included positions in core genome","Percentage of valid and included positions in core genome", "Percentage of all positions that are valid, included, and part of the core genome");

    print $out_fh join("\t",@header) . "\n";
    #negative final_invalid indicates there was no invalid positions file given
    my ($final_core,$final_total,$final_invalid)= (0,0,0);

    if (not  $invalid_total) 
    {
        $final_invalid='N/A';
    }

    foreach my $chrom( sort {$a cmp $b } keys %vcfcore) 
    {
        my ($core,$total) = ($vcfcore{$chrom}{'core'},$vcfcore{$chrom}{'total'});
        my $invalid = 'N/A';
        my $total_no_invalid;

        my $perc;

        #if there was invalid found
        if (exists $invalid_total->{$chrom}) 
        {
            $invalid = $invalid_total->{$chrom}; #grab the total number of invalid positions marked by the user for the current chrom
            $final_invalid +=$invalid;

            if ($total == $invalid ) 
            {
                $perc ="NaN";
                $total_no_invalid=0;
            }
            else 
            {
                $perc = sprintf("%.2f",($core/( $total-$invalid)*100));
                $total_no_invalid = $total-$invalid;
            }
        }
        elsif ( $vcfcore{$chrom}{'invalid'} && not exists $invalid_total->{$chrom}) 
        {
            die "Found invalid positions for '$chrom' but was not found in invalid positions file.\n";
        }
        else 
        {
            $perc = sprintf("%.2f",($core/$total*100));
            $total_no_invalid = $total;
        }
        $final_core +=$core;
        $final_total +=$total;
        my $final_totals_with_invalid_contigs = sprintf("%.2f",($core/$total)*100);
        print $out_fh join ("\t", ($chrom,$total,$invalid,$total_no_invalid,$core,$perc,$final_totals_with_invalid_contigs)) . "\n";
    }


    #getting the total for all references
    my ($perc,$final_total_no_invalid);

    if ($final_invalid eq 'N/A' || $final_invalid == 0 ) 
    {
        $perc = sprintf("%.2f",$final_core/$final_total*100);
        $final_total_no_invalid=$final_total;
    }
    else 
    {
        if ($final_total == $final_invalid ) 
        {
            $perc ="NaN";
            $final_total_no_invalid=0;
        }
        else 
        {
            $perc = sprintf("%.2f",$final_core/($final_total - $final_invalid)*100);
            $final_total_no_invalid=$final_total-$final_invalid;
        }
    }
    my $final_totals_with_invalid = sprintf("%.2f",($final_core/$final_total)*100);
    print $out_fh join ("\t",'all',$final_total,$final_invalid,$final_total_no_invalid,$final_core,$perc,$final_totals_with_invalid) . "\n";

    return;
}


    sub prepare_inputs 
{
    my (%consolidate_vcf, @formats, $output_base, $reference, $fasta, $invalid, $requested_cpus, $bcftools);
		my ($refs_info, $invalid_pos, $invalid_total, $help, $verbose, $inc_list, $exc_list, $formats);

    if( @_ && $_[0] eq __PACKAGE__ )
	{
		GetOptions(
		     'consolidate_vcf=s'     => \%consolidate_vcf,
				 'format|f=s'            => \@formats,
				 'output-base|o=s'       => \$output_base,
				 'reference|r=s'         => \$reference,
				 'fasta=s'               => \$fasta,
				 'invalid-pos=s'         => \$invalid,
				 'numcpus=i'             => \$requested_cpus,
				 'b|bcftools-path=s'     => \$bcftools,
				 'help|h'                => \$help,
				 'verbose|v'             => \$verbose,
				 'include=s'             => \$inc_list,
 				 'exclude=s'             => \$exc_list
		);
		pod2usage(0) if $help;
	}
	else
	{
		my $vcfs;
	    ($vcfs,$formats,$output_base,$reference,$fasta,$invalid,$requested_cpus,$bcftools,$inc_list,$exc_list) = @_;
		%consolidate_vcf = %{$vcfs};
	}

    $verbose = 0 if (not defined $verbose);

    if ( scalar keys %consolidate_vcf == 0)
    {
        print STDERR "Was not able to find any vcf files from consolidate_vcf.\n\n";
		pod2usage(1);
    }

	if(not defined $output_base)
	{
		print STDERR "output-base undefined\n";
		pod2usage(1);
	}

	if (not defined $bcftools)
	{
        #check to see if bcftools is on the path
        #normally we always want to be passed the path to the tool but this is for Galaxy implementation
        my $alive=`bcftools 2>&1 1>/dev/null`;
        if ( $alive && $alive =~ /Program: bcftools/) 
        {
            $bcftools="bcftools";
        }
        else 
        {
            print STDERR "bcftools-path not defined and not found on the path.\n";
			pod2usage(1);
        }
    }

    #need check to see if bcftools was complied with htslib and also has the correct plugin installed
	
	my $usage_state = `$bcftools 2>&1 1>/dev/null`;
    if ( not $usage_state =~ /Version: .* \(using htslib/ ) 
    {
        print STDERR "bctools was not complied with htslib.\nPlease re-compile with htslib\nInstruction: http://samtools.github.io/bcftools/\n";
		pod2usage(1);
    }

    if (not defined $reference)
    {
        print STDERR "reference name not defined, calling it 'reference'\n";
        $reference = 'reference';
    }

    if (defined $invalid)
    {
        if ( ! -e $invalid){
            print STDERR "Was given an invalid position file but could not locate it '$invalid'\n";
			pod2usage(1);
        }
    }
    else
    {
        print "invalid position file not defined, Will ignore step\n";
    }

	$requested_cpus = 1 if (not defined $requested_cpus);

    if (@formats <= 0)
    {
		print STDERR "warning: format not defined, assuming fasta\n";
		@formats = ("fasta");
    }

    else
    {
        for my $format (@formats)
        {
			if(not defined $valid_formats{$format})
			{
				print STDERR "unrecognized format '$format', must be one of '".join(' ', keys %valid_formats),"'\n";
				pod2usage(1);
			}
        }
    }

    if ( not -e $fasta) 
    {
        print STDERR "Error: Was not given reference fasta file\n";
		pod2usage(1);
    }

  	$refs_info = refs_info($fasta);

    if ($invalid)
    {
        my $invalid_positions_parser = InvalidPositions->new;
        ($invalid_pos,$invalid_total) = $invalid_positions_parser->read_invalid_positions($invalid);
    }

	if ($inc_list && $exc_list) 
	{
        die "Cannot have both an include and exclude list. Please specify only one!\n";
    }

    #if we have a including list, then use it
    if ($inc_list) 
    {
        %consolidate_vcf = %{strain_selection($inc_list,1,\%consolidate_vcf)};
    }

	#if we have a excluding list, then use it
    if ($exc_list) 
    {
        %consolidate_vcf = %{strain_selection($exc_list,0,\%consolidate_vcf)};
    }

	return (\%consolidate_vcf,\@formats,$output_base,$reference,$invalid_pos,$invalid_total,$requested_cpus,$bcftools,$refs_info);
}

sub strain_selection 
{
	my ($list,$keep,$consolidate_vcf)=@_;

	my $vcfs = $consolidate_vcf;
	my %tokeep_vcf;

	#get list of all strains and make it lower case
	#only need to check one of the two hashes because both should be the same by this point
	my %strains = map { lc $_ => $_ } keys %{$consolidate_vcf};


	open my $fh, '<',$list || die "Could not open file '$list'\n";
	while (my $name = <$fh>) 
	{
			chomp $name;
			$name = lc $name; #lc so we can do case insensitive lookup

			#check to see if the name is in the list of keys
			if ( exists $strains{$name}) 
			{
					#if we are doing a inclusion list, add it to the new hashes
					#otherwise just delete from the hash coming in
					my $org_name = $strains{$name};
					if ($keep) 
					{
						$tokeep_vcf{$org_name} = $consolidate_vcf->{$org_name};
					}
					else 
					{
						delete $vcfs->{$org_name};
					}
			}
			else 
			{
				print STDERR "WARNING: Could not find strain '$name' in list of strains. Ignorning\n";
			}
	}
	close $fh;

	#if we are doing an exclusion list, then move return hashes
	if (not $keep) 
	{
		return $vcfs;
	}
	else 
	{
		return \%tokeep_vcf;
	}
}

1;

=pod

=head1 NAME

quick2snv_alignment.pl

=head1 SYNOPSIS

quick2snv_alignment.pl --consolidate_vcf v1=files/dataset1.dat --consolidate_vcf v2=files/dataset2.dat --consolidate_vcf v3=files/dataset3.dat --format fasta --format phylip --output-base /tmp/results --reference strain_24 --fasta /files/reference.fasta --invalid-pos [invalid positions TSV file] --numcpus 5  --bcftools-path /opt/bcftools/bcftools

=head1 OPTIONS

=over

=item B<--consolidate_vcf> [REQUIRED]

Hash containing combined vcf files from consolidate_vcfs.

=item B<--format> [OPTIONAL]

The format to output the alignment to, one of the Bio::AlignIO supported formats (default: fasta).

=item B<--output-base> [REQUIRED]

The output base name for the alignment file(s).

=item B<--reference> [OPTIONAL]

The name of the reference to use in the alignment (default: reference).

=item B<--fasta> [REQUIRED]

Fasta file.

=item B<--invalid-pos> [OPTIONAL]

A TSV file that contains a list of range(s) (one per line) of CHROM\\tSTART_POS\\tEND_POS\\n".

=item B<--numcpus> [REQUIRED]

Desired number of CPUs for the job.

=item B<--bcftools-path> [OPTIONAL]

Path to BCFTools.

=item B<-h>, B<--help>

Displays the help screen.

=back

=head1 DESCRIPTION

Given a hash of consolidated vcfs, examines all called SNPs and generates a 'snv_alignment' of all the high quality SNPs/base calls (to be used for further phylogenetic analysis).

=cut
