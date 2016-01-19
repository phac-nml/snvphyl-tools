#!/usr/bin/env perl

# vcf2snv_alignment
# Purpose:  Given a set of combined freebayes and mpileup files, examines all called SNPs and generates a 'snv_alignment'
#  of all the high quality SNPs/base calls (to be used for further phylogenetic analysis).

use warnings;
use strict;
use FindBin;
use lib $FindBin::Bin.'/lib';
use lib $FindBin::Bin;
use Pod::Usage;
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


	#create alignment files
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

sub filter_positions 
{
    my ($consolidate_vcf,$refs,$invalid_pos,$valid_positions,$bcftools,$cpus) = @_;

    my %results;
    my $pm;
    my $num_cpus=`cat /proc/cpuinfo | grep processor | wc -l`;
    chomp $num_cpus;
    #ensure that you user cannot request more threads then CPU on the machine
    if ( $cpus > $num_cpus)
    {
        $cpus = $num_cpus;
    }

    my %vcfcore;


    $pm=Parallel::ForkManager->new($cpus);

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

    my $bit_size = 100000;
    my $job_id=0;
    my $tmp_dir = tempdir (CLEANUP => 1);

    #print header file
    open my $out, '>', "$tmp_dir/$job_id";
    print $out "#Chromosome\tPosition\tStatus\tReference\t";
    my @samples_list = sort {$a cmp $b } keys %{$consolidate_vcf};
    print $out join("\t",@samples_list);
    print $out "\n";
    close $out;


    foreach my $chrom( sort { $a cmp $b } keys %$refs) 
    {
        my ($start,$stop)=(0,0);
        my ($range);
        my $length = $refs->{$chrom}{'length'};
        $vcfcore{$chrom}= {'invalid' => 0,
                           'core' => 0,
                           'total'=>$length
                          };

        if ( !$length) 
        {
            die "Could not determine length of '$chrom' from provided reference fasta file\n";
        }


        while ($stop < $length ) 
        {
            #inclusive range, we do not care if we go over b/c vcf query will only return what it can.
            $start = $stop +1;
            $stop = $stop+$bit_size;
            $range="$chrom:" . join('-',($start,$stop));

            $range = quotemeta $range;

            $job_id++;

            my $pid = $pm->start and next if $parallel;

            my %stats = ($chrom => {'invalid' => 0,
                                    'core' => 0,
                                   });

            my $f_name = "$tmp_dir/$job_id";
            open my $out, '>',$f_name;

            my $streamers = Streaming::create_streamers($consolidate_vcf,$range,$job_id,$bcftools);
            my $cur_pos = $start;

            my @data = $streamers->($chrom,$cur_pos);


            #search for all entries to have "EOF" as their data.
            while (any { $_->{'status'} ne 'EOF' } @data) 
            {
                #get current bp for current positions
                my $ref_bp = $refs->{$chrom}{'bps'}[$cur_pos];

                #if we do not have any SNP, we do not print at all!
                if ( any {  exists $_->{'alt'} && $_->{'alt'} ne '.' } @data ) 
                {
                    my @line = ($chrom,$cur_pos);

                    #all columns are passed all cut off parameters
                    if ( all { $_->{'status'} eq 'PASS' }  @data) 
                    {
                        if ($invalid_pos && exists $invalid_pos->{"${chrom}_${cur_pos}"} ) 
                        {
                            push @line, 'filtered-invalid';
                            $stats{$chrom}{'invalid'}++;
                        }
                        else 
                        {
                            push @line, 'valid';
                            $stats{$chrom}{'core'}++;
                        }

                        push @line,$ref_bp;
                        foreach ( @data) 
                        {
                            if ( $_->{'alt'} eq '.') 
                            {
                                push @line,$ref_bp;
                            }
                            else 
                            {
                                push @line,$_->{'alt'};
                            }
                        }
                        print $out join("\t",@line) . "\n";

                    }
                    #now we are working with position where there is at least one SNP and at least one filtered-*
                    elsif (any {  ( $_->{'status'} eq 'PASS' || $_->{'status'} eq 'filtered-invalid' )  && $_->{'alt'} ne '.'} @data )
                    {
                        my $is_core=1; #if all positions have coverage, part of the core regardless if we do not agree it's a high quality SNP or invalid position
                        if ($invalid_pos && exists $invalid_pos->{"${chrom}_${cur_pos}"} ) 
                        {
                            push @line, 'filtered-invalid';
                            $stats{$chrom}{'invalid'}++;
                            $is_core=0;
                        }
                        else 
                        {
                            #always take filtered-invalid before anything else if it there.
                            my $index = firstidx {$_->{'status'} eq 'filtered-invalid'  } @data;

                            #if we did not find a filtered-invalid, find another none valid status
                            if ( $index == -1) 
                            {
                                $index = firstidx {$_->{'status'} ne 'PASS' } @data;
                            }

                            if ( $data[$index]->{'status'} eq '' || $data[$index]->{'status'} eq 'EOF') 
                            {
                                push @line,'filtered-coverage';
                            }
                            else 
                            {
                                my $stats =  $data[$index]->{'status'};
                                $stats = 'filtered-' . $stats;
                                push @line,$stats;
                            }
                        }

                        push @line,$ref_bp;

                        foreach my $col ( @data) 
                        {
                            my $status = $col->{'status'};

                            #if '', implies there is no reads covering that $cur_pos
                            if ( $status eq '' || $status eq 'EOF' || $status eq 'coverage') 
                            {
                                push @line,'-';
                                $is_core=0;
                            }
                            else 
                            {
                                #if we have filtered-mpileup, we have inconsistent calles between variant callers
                                if ($status eq 'mpileup' ) 
                                {
                                    push @line,'N';
                                }
                                #have a position that passes the cut-off parameter. Either show the SNP or the reference
                                elsif ( $status eq 'PASS' || $status eq 'filtered-invalid') 
                                {
                                    $is_core=0 if $status eq 'filtered-invalid';

                                    if ( $col->{'alt'} eq '.') 
                                    {
                                        push @line,$ref_bp;
                                    }
                                    else 
                                    {
                                        push @line,$col->{'alt'};
                                    }
                                }
                                else 
                                {
                                    push @line,'-';
                                }

                            }
                        }

                        $stats{$chrom}{'core'}++ if $is_core; #if all positions had coverage, add to the core

                        print $out join("\t",@line) . "\n";

                    }#end else
                    else 
                    {
                        print "WARNING: Edge case found at '$chrom' on position '$cur_pos'\n" if $verbose;
                    }

                }#end if
                
                #see if everyone has a PASS status and there is no SNP
                elsif ( all { $_->{'status'} eq 'PASS' }  @data ) 
                {
                    #check to see if the positions in the invalid position, if not, count it toward the total
                    if (not ($invalid_pos && exists $invalid_pos->{"${chrom}_${cur_pos}"} )) 
                    {
                        $stats{$chrom}{'core'}++;
                    }


                }
                else 
                {
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

    my @header = ("#Reference","total length","total invalid pos","total valid pos","total core","Percentage in core");

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

        print $out_fh join ("\t", ($chrom,$total,$invalid,$total_no_invalid,$core,$perc)) . "\n";
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
    print $out_fh join ("\t",'all',$final_total,$final_invalid,$final_total_no_invalid,$final_core,$perc) . "\n";

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
		pod2usage(1) if $help;
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

vcf2snv_alignment.pl

=head1 SYNOPSIS

vcf2snv_alignment.pl --consolidate_vcf v1=files/dataset1.dat --consolidate_vcf v2=files/dataset2.dat --consolidate_vcf v3=files/dataset3.dat --format fasta --format phylip --output-base /tmp/results --reference strain_24 --fasta /files/reference.fasta --invalid-pos [invalid positions TSV file] --numcpus 5 --bcftools-path /opt/bcftools/bcftools

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
