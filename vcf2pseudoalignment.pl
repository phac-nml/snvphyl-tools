#!/usr/bin/env perl
# vcf2pseudoalignmnet
# Purpose:  Given a set of *.vcf files, examines all called SNPs and generates a 'pseudoalignment'
#  of all the high quality SNPs/base calls (to be used for further phylogenetic analysis).

use warnings;
use strict;

use Getopt::Long;
use Storable qw /dclone store retrieve/;
use File::Basename;
use Parallel::ForkManager;

use Bio::AlignIO;
use Bio::SimpleAlign;
# VCF module from vcftools: http://vcftools.sourceforge.net/index.html
use Vcf;
use File::Temp qw /tempdir/;


my $verbose;
my $unknown_base = 'N';

my $snp_info = {'removed' => {'insertions' => 0, 'deletions' => 0, 'multi' => 0, 'other' => 0},
		'total' => 0,
		'positions' => 0,
		'snps' => {'filtered-coverage' => 0, 'filtered-mpileup' => 0, 'kept' => 0}
		};

sub usage
{
	"Usage: $0 --vcf-dir [vcf dir] --mpileup-dir [mpileup dir] --output-base [base of output alignment file]\n".
	"Parameters:\n".
	"\t--vcf-dir: The directory containing the vcf files.\n".
	"\t--mpileup-dir: Directory containing the vcf files produced by 'samtools mpileup \$file | bcftools view -cg' (*.vcf).\n".
	"\t-o|--output-base:  The output base name for the alignment file(s)\n".
	"Options:\n".
	"\t-u|--uniquify:  Make the seq names unique and print mapping file to real names (for phylip format limitations)\n".
	"\t-r|--reference:  The name of the reference to use in the alignment (default: reference)\n".
	"\t-f|--format:  The format to output the alignment to, one of the Bio::AlignIO supported formats (default: fasta)\n".
	"\t-c|--coverage-cutoff:  The cutoff for coverage to include a reference base (default: 1)\n".
	"\t--invalid-pos: A TSV file that contains a list of range(s) (one per line) of CHROM\\tSTART_POS\\tEND_POS\\n".
	"\t--verbose:  More information printed\n".
	"\t--keep-ambiguous:  Keep ambiguous characters in alignment output file\n".
	"\t-h|--help:  Help\n";
}



sub parse_invalid 
{
    my ($file) =  @_;
    my %invalid;

    open(my $fh, "<" , "$file") or die "Could not open $file: $!";

    while(my $line = readline($fh))
    {
	chomp $line;
	my ($sub_line) = ($line =~ /^([^#]*)/);
	my ($chrom,$start,$end) = split(/\t/,$sub_line);
	next if (not defined $chrom or $chrom eq '');
	next if ($start !~ /^\d+$/);
	next if ($end !~ /^\d+$/);

	# swap in case start/end are reversed
	my $real_start = ($start < $end) ? $start : $end;
	my $real_end = ($start < $end) ? $end : $start;


        foreach my $i ( $real_start..$real_end ) {
	    $invalid{"${chrom}_${i}"} = 1;
        }
    }

    close($fh);
    return  \%invalid
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

# Parses position data structure (for a single chromosome) and builds alignment.
# Input:  positions_hash  The hash of all positions of variants
#	  samples_list  A list of all sample names to proccess
# Output:  An alignment of all valid variant positions (in FASTA format)
sub variants_alignment
{
	my ($positions_hash, $chromosome, $reference, $samples_list, $mpileup_data, $coverage_cutoff,$invalid_pos,$keep_ambiguous) = @_;
	my $alignment_string = undef;

	# stores pseudoalignment in form of
	# {sample => {alignment => 'string', positions => [pos array]}}
	my %alignment;

	# for printing out list of valid/excluded positions
	# form of {pos => {status => valid/invalid, ref => 'base', 'samples' => {sample => 'base'}}}
	my $total_positions = {};
	for my $sample (@$samples_list)
	{
		$alignment{$sample} = {'alignment' => '', 'positions' => []};
	}

	# add alignment string for reference
	$alignment{$reference} = {'alignment' => '', 'positions' => []};

	for my $pos (sort {$a <=> $b} keys %$positions_hash)
	{
		my $sample_hash = $positions_hash->{$pos};
		my @sample_hash_list = keys %$sample_hash;
		my $first_sample = $sample_hash->{$sample_hash_list[0]};
		my $ref_base = $first_sample->{'ref'}; # get reference base
		$total_positions->{$pos} = {'ref' => $ref_base} if (not defined $total_positions->{$pos});

		# check for which samples have no variant called in this position
		my $alignment_local = {};
		for my $sample (@$samples_list)
		{
			# get data for aligment we are building
			if (not exists $sample_hash->{$sample})
			{
				# check for samtools mpileup of this position
				my $pileup_vcf = $mpileup_data->{$sample}->{$chromosome}->{$pos};
				if (not defined $pileup_vcf)
				{
					print STDERR "fail for $sample:$chromosome:$pos, mpileup(coverage) not defined\n" if ($verbose);
					$alignment_local->{$sample} = {'base' => $unknown_base, 'position' => $pos};
					$snp_info->{'snps'}->{'filtered-coverage'}++;

					$total_positions->{$pos}->{'samples'}->{$sample} = '-';
					my $is_valid = $total_positions->{$pos}->{'status'};
					if (not defined $is_valid or $is_valid eq 'valid')
					{
						$total_positions->{$pos}->{'status'} = 'filtered-coverage';
					}
				}
				else
				{
					my $ref = $pileup_vcf->{'ref'};
					my $alt = $pileup_vcf->{'alt'};
					my $coverage = $pileup_vcf->{'cov'};

					if (not defined $coverage)
					{
						die "Error: mpileup coverage not defined for sample $sample:$chromosome:$pos\n";
					}
					elsif (not defined $alt)
					{
						die "Error: mpileup for $sample:$chromosome:$pos, alt not defined\n" if ($verbose);
					}
					elsif (not defined $ref)
					{
						die "Error: mpileup for $sample:$chromosome:$pos, ref not defined\n" if ($verbose);
					}
					elsif ($coverage < $coverage_cutoff)
					{
						print STDERR "fail for $sample:$chromosome:$pos, coverage=$coverage < cutoff=$coverage_cutoff\n" if ($verbose);
						$alignment_local->{$sample} = {'base' => $unknown_base, 'position' => $pos};
						$snp_info->{'snps'}->{'filtered-coverage'}++;

						$total_positions->{$pos}->{'samples'}->{$sample} = '-';
						my $is_valid = $total_positions->{$pos}->{'status'};
						if (not defined $is_valid or $is_valid eq 'valid')
						{
							$total_positions->{$pos}->{'status'} = 'filtered-coverage';
						}
					}
					elsif ($alt ne '.')
					{
						print STDERR "fail for $sample:$chromosome:$pos, mpileup data gives alt=$alt (ref=$ref), but no variant called from vcf data\n" if ($verbose);
						$alignment_local->{$sample} = {'base' => $unknown_base, 'position' => $pos};
						$snp_info->{'snps'}->{'filtered-mpileup'}++;

						$total_positions->{$pos}->{'samples'}->{$sample} = $unknown_base;
						my $is_valid = $total_positions->{$pos}->{'status'};
						if (not defined $is_valid or $is_valid eq 'valid')
						{
							$total_positions->{$pos}->{'status'} = 'filtered-mpileup';
						}
					}
					elsif ($invalid_pos && exists $invalid_pos->{"${chromosome}_${pos}"} ) 
					{
					    print STDERR "fail for $sample:$chromosome:$pos, position is invalid due to invalid position file provided by operator\n" if ($verbose);
					    $snp_info->{'snps'}->{'filtered-invalid'}++;

					    $total_positions->{$pos}->{'samples'}->{$sample} = $ref_base;
					    $total_positions->{$pos}->{'status'} = 'filtered-invalid';
					}
					else
					{
						print STDERR "pass for $sample:$chromosome:$pos, coverage=$coverage > cutoff=$coverage_cutoff and no variant called from mpileup data (alt=$alt, ref=$ref)\n" if ($verbose);
						$alignment_local->{$sample} = {'base' => $ref_base, 'position' => $pos};

						$snp_info->{'snps'}->{'kept'}++;
						$total_positions->{$pos}->{'samples'}->{$sample} = $ref_base;
						my $is_valid = $total_positions->{$pos}->{'status'};
						if (not defined $is_valid)
						{
							$total_positions->{$pos}->{'status'} = 'valid';
						}
					}
				}
			}
			else
			{
				my $pileup_vcf = $mpileup_data->{$sample}->{$chromosome}->{$pos};
				die "Error: mpileup for $sample:$chromosome:$pos not defined, but position called in other alignment software" if (not defined $pileup_vcf);
				my $coverage = $pileup_vcf->{'cov'};
				die "Error: position defined, but coverage in mpileup for $sample:$chromosome:$pos not defined\n" if (not defined $coverage or $coverage !~ /^\d+$/);
				if ($coverage < $coverage_cutoff)
				{
					print STDERR "fail for $sample:$chromosome:$pos, coverage=$coverage < cutoff=$coverage_cutoff\n" if ($verbose);
					$alignment_local->{$sample} = {'base' => $unknown_base, 'position' => $pos};
					$snp_info->{'snps'}->{'filtered-coverage'}++;

					$total_positions->{$pos}->{'samples'}->{$sample} = '-';
					my $is_valid = $total_positions->{$pos}->{'status'};
					if (not defined $is_valid or $is_valid eq 'valid')
					{
						$total_positions->{$pos}->{'status'} = 'filtered-coverage';
					}
				}
				elsif ($invalid_pos && exists $invalid_pos->{"${chromosome}_${pos}"} ) 
				{
				    print STDERR "fail for $sample:$chromosome:$pos, position is invalid due to invalid position file provided by operator\n" if ($verbose);
				    $snp_info->{'snps'}->{'filtered-invalid'}++;
				    
				    $total_positions->{$pos}->{'samples'}->{$sample} = $sample_hash->{$sample}->{'alt'};
				    $total_positions->{$pos}->{'status'} = 'filtered-invalid';
				}
				else
				{
					$alignment_local->{$sample} = {'base' => $sample_hash->{$sample}->{'alt'}, 'position' => $pos};
					$snp_info->{'snps'}->{'kept'}++;
					$total_positions->{$pos}->{'samples'}->{$sample} = $sample_hash->{$sample}->{'alt'};
					my $is_valid = $total_positions->{$pos}->{'status'};
					if (not defined $is_valid)
					{
						$total_positions->{$pos}->{'status'} = 'valid';
					}
				}
			}
		}

		# fill in overall data structure
		if (($keep_ambiguous or $total_positions->{$pos}->{'status'} eq 'valid')
			 and $total_positions->{$pos}->{'status'} ne 'filtered-invalid')
		{
			for my $sample (keys %$alignment_local)
			{
				my $base = $alignment_local->{$sample}->{'base'};
				my $pos = $alignment_local->{$sample}->{'position'};

				my $alignment_sample = $alignment{$sample};

				$alignment_sample->{'alignment'} .= $base;
				push(@{$alignment_sample->{'positions'}}, $pos);
			}

			# fill in data for reference
			$alignment{$reference}->{'alignment'} .= $ref_base;
			push(@{$alignment{$reference}->{'positions'}}, $pos);
		}
		else
		{
			print STDERR "skipping over position $chromosome:$pos with status ".$total_positions->{$pos}->{'status'}."\n" if ($verbose);
		}
	}


	return (\%alignment,$total_positions);
}

# parse_variants
# Purpose: parses variant files and fills in data stucture.
# data structure stores data in format given below
# 'chrom' => { position => { sample_name => 
#			{ ref => ref_base, alt => alt_base}
#		   }
#	     }
#
# Input: vcf_files  A hash giving the vcf files to parse
# Output:  A reference to the vcf_data structure
sub parse_variants
{
	my ($vcf_files,$requested_cpus) = @_;
        
         my $pm;
         my $num_cpus=`cat /proc/cpuinfo | grep processor | wc -l`;
         chomp $num_cpus;
         #ensure that you user cannot request more threads then CPU on the machine
         if ( $requested_cpus > $num_cpus) {
             $requested_cpus = $num_cpus;
         }
         my %list;
        
         $pm=Parallel::ForkManager->new($requested_cpus);
         # data structure retrieval and handling
          $pm -> run_on_finish ( # called BEFORE the first call to start()
              sub {
                  my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $child_data) = @_;
                  # retrieve data structure from child
                  if (defined($child_data)) {  # children are forced to send anything
                      my ($single_key) = keys %{$child_data};
                      $list{$single_key} = $child_data->{$single_key};
                  } else {
                      die "One or more mpileup file did not produce any data!\n";
                  }
              }
          );

        #creating temporary directory on local machine to store results from each child.
        my $tmpdir = tempdir( CLEANUP => 1 );

	my %vcf_data;
	for my $vcf_key (keys %$vcf_files)
	{
                $pm->start and next;
                my %vcf_data;
		my $vcf_file = $vcf_files->{$vcf_key};

		print STDERR "Working on $vcf_file\n" if ($verbose);
		my $vcf = Vcf->new(file => $vcf_file);
		$vcf->parse_header();
		while (my $data = $vcf->next_data_array)
		{
			if (@$data <= 4)
			{
				die "Not enough information for vcf file $vcf_file, line ".join(' ',@$data);
			}
			else
			{
				my $chrom = $data->[0];
				my $position = $data->[1];
				my $ref = $data->[3];
				my $alt = $data->[4];
	
				$snp_info->{'total'}++;
				# check for indel/multi-snp
				# case: multi-snp
				if ((length($ref) > 1) and length($ref) eq length($alt))
				{
                                        #print STDERR "multi variant region found in $vcf_file, skipping: \"".join(' ',@$data),"\"\n" if ($verbose);
					$snp_info->{'removed'}->{'multi'}++;
				}
				elsif (length($ref) > 1 and length($ref) < length($alt))
				{
					#print STDERR "length($ref) > 1 and length($ref) < length($alt) in $vcf_file, skipping: \"".join(' ',@$data),"\"\n" if ($verbose);
					$snp_info->{'removed'}->{'other'}++;
				}
				elsif (length($ref) > 1) # deletion in query strain
				{
					$snp_info->{'removed'}->{'deletions'}++;
					#print STDERR "deletion found in $vcf_file, skipping: \"".join(' ',@$data),"\"\n" if ($verbose);
				}
				elsif (length($alt) > 1) # insertion in query strain
				{
					$snp_info->{'removed'}->{'insertions'}++;
					#print STDERR "insertion found in $vcf_file, skipping: \"".join(' ',@$data),"\"\n" if ($verbose);
				}
				elsif((length($alt) == 1) and length($ref) == 1) # SNP
				{
					#print STDERR "for file $vcf_file, keeping variant \"".join(' ',@$data),"\"\n" if ($verbose);
					my $chrom_hash;
					if (not exists $vcf_data{$chrom})
					{
						$chrom_hash = {};
						$vcf_data{$chrom} = $chrom_hash;
					}
					else
					{
						$chrom_hash = $vcf_data{$chrom};
					}
		
					my $position_hash;
					if (not exists $chrom_hash->{$position})
					{
						$position_hash = {};
						$chrom_hash->{$position} = $position_hash;
						$snp_info->{'positions'}++;
					}
					else
					{
						$position_hash = $chrom_hash->{$position};
					}
		
					die "Error: duplicate position $position for vcf file $vcf_file" if (exists $position_hash->{$vcf_key});
					$position_hash->{$vcf_key} = {'ref' => $ref, 'alt' => $alt};
				}
				else # any other case?
				{
					#print STDERR "invalid lengths for ref and query found in $vcf_file, skipping: \"".join(' ',@$data),"\"\n" if ($verbose);
					$snp_info->{'removed'}->{'other'}++;
				}

			}
                    }
                

                #storing
                #should be a single key only since we having one iteration per thread
                store \%vcf_data, "$tmpdir/$vcf_key";
                $pm->finish(0,{$vcf_key=>"$tmpdir/$vcf_key"});
            }
        
        $pm->wait_all_children;
        
        #read all the files from dir storable into %vcf_data and go on!
        foreach my $vcf_key( keys %list) {
            my $data = retrieve($list{$vcf_key});
            #go thru each chromsome and add the sample
            foreach my $chrom (keys %$data ) {
                foreach my $position(keys %{$data->{$chrom}} ) {
                    $vcf_data{$chrom}{$position}{$vcf_key} = $data->{$chrom}{$position}{$vcf_key};
                }

            }

        }
        

	return \%vcf_data;
}

sub parse_mpileup
{
	my ($mpileup_files, $vcf_data,$requested_cpus) = @_;
	# vcf_data is hash constructed from parse_variants
	# converts to a positions_hash mapping 'chrom' => {'pos' => undef} for all positions
	my $positions_hash = {};
	for my $chr (keys %$vcf_data)
	{
		my $chr_hash = $vcf_data->{$chr};
		my $new_chr_hash = {};
		for my $pos (keys %$chr_hash)
		{
			$new_chr_hash->{$pos} = undef;
		}
		$positions_hash->{$chr} = $new_chr_hash;
	}
	# positions_hash gets cloned and filled in for every sample
	# fills in to 'chrom' => {'pos' => vcf-info/undef}

	my %sample_pileup_hash;
        my %list;
        
        my $pm;
        my $num_cpus=`cat /proc/cpuinfo | grep processor | wc -l`;
        chomp $num_cpus;
        #ensure that you user cannot request more threads then CPU on the machine
        if ( $requested_cpus > $num_cpus) {
            $requested_cpus = $num_cpus;
            
        }
        $pm=Parallel::ForkManager->new($requested_cpus);
        # data structure retrieval and handling
         $pm -> run_on_finish ( # called BEFORE the first call to start()
             sub {
                 my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $child_data) = @_;
                 # retrieve data structure from child
                 if (defined($child_data)) {  # children are forced to send anything
                     my ($single_key) = keys %{$child_data};
                     $list{$single_key} = $child_data->{$single_key};
                 } else {
                     die "One or more mpileup file did not produce any data!\n";
                 }
             }
         );

        #creating temporary directory on local machine to store results from each child.
        my $tmpdir = tempdir( CLEANUP => 1 );
        
	# read through each vcf file
	# place into mpileup_hash, keeping only vcf lines necessary (in sorted_positions)

	for my $sample (keys %$mpileup_files)
	{
                $pm->start and next;
		my $working_mpileup = dclone($positions_hash);
		my $vcf_file = $mpileup_files->{$sample};

		print STDERR "reading mpileup information from $vcf_file\n" if ($verbose);
		my $vcf = Vcf->new(file => $vcf_file) or die "Could not parse vcf file $vcf_file";
		$vcf->parse_header();
		while (my $data = $vcf->next_data_array)
		{
			if (@$data <= 4)
			{
				die "Not enough information for vcf file $vcf_file, line ".join(' ',@$data);
			}
			else
			{
				my $chrom = $data->[0];
				my $position = $data->[1];
				
				# found a position to keep
				if (exists $working_mpileup->{$chrom}->{$position})
				{
					my %info_hash = %{variant_info_to_hash($data->[7])};
					my $ref = $data->[3];
					my $alt = $data->[4];
					my $coverage = $info_hash{'DP'};

					if (exists $info_hash{'INDEL'})
					{
						#print STDERR "skipping INDEL vcf-line found for $sample: '".join(" ",@$data)."\n";
					}
					elsif (defined $working_mpileup->{$chrom}->{$position})
					{
						die "Duplicate position found in $sample.\nvcf-line: ".join(" ",@$data);
					}
					else
					{
                                                #print STDERR "mpileup info found for $sample:$chrom:$position\n" if ($verbose);
						$working_mpileup->{$chrom}->{$position} = {'ref' => $ref, 'alt' => $alt, 'cov' => $coverage};
					}
				}
			}
                    }
            
                $vcf->close();

                #storing
                my $name = basename($sample);
                store $working_mpileup, "$tmpdir/$name";
                
                $pm->finish(0,{$sample=>"$tmpdir/$name"})

            
	}
        
        $pm->wait_all_children;

         #read all the files from dir storable into sample_pileup_hash and go on!
         foreach my $sample( keys %list) {
             $sample_pileup_hash{$sample} = retrieve($list{$sample});
         }
	return \%sample_pileup_hash;
}

############
### MAIN ###
############

# maps format name to format file extension
my %valid_formats = ('fasta' => 'fasta', 'phylip' => 'phy', 'clustalw' => 'cl');

my $vcf_dir;
my $mpileup_dir;
my $output_base;
my @formats;
my $reference;
my $coverage_cutoff;
my $uniquify;
my $help;
my $requested_cpus;
my $invalid;
my $keep_ambiguous;

my $command_line = join(' ',@ARGV);

if (!GetOptions('vcf-dir|d=s' => \$vcf_dir,
		'mpileup-dir|b=s' => \$mpileup_dir,
		'format|f=s' => \@formats,
		'output-base|o=s' => \$output_base,
		'reference|r=s' => \$reference,
		'coverage-cutoff|c=i' => \$coverage_cutoff,
		'uniquify|u' => \$uniquify,
		'invalid-pos=s' => \$invalid,
		'help|h' => \$help,
		'keep-ambiguous' => \$keep_ambiguous,
                'numcpus=i' => \$requested_cpus,
                
		'verbose|v' => \$verbose))
{
	die "Invalid option\n".usage;
}

print usage and exit(0) if (defined $help);
$verbose = 0 if (not defined $verbose);

die "vcf-dir undefined\n".usage if (not defined $vcf_dir);
die "vcf-dir does not exist\n".usage if (not -e $vcf_dir);

die "mpileup-dir undefined\n".usage if (not defined $mpileup_dir);
die "mpileup-dir does not exist\n".usage if (not -e $mpileup_dir);

die "output-base undefined\n".usage if (not defined $output_base);

if (not defined $reference)
{
	print STDERR "reference name not defined, calling it 'reference'\n";
	$reference = 'reference';
}



if (defined $invalid)
{
    if ( ! -e $invalid)
    {
	die "Was given an invalid position file but could not locate it '$invalid'\n";
    }
}
else
{
    print STDERR "invalid position file not defined, Will ignore step\n";
}


$requested_cpus = 1 if (not defined $requested_cpus);
$uniquify = 0 if (not defined $uniquify);

if (@formats <= 0)
{
	print STDERR "warning: format not defined, assuming fasta\n";
	@formats = ("fasta");
}
else
{
	for my $format (@formats)
	{
		die "unrecognized format '$format', must be one of '".join(' ', keys %valid_formats),"'\n" if (not defined $valid_formats{$format});
	}
}

if (not defined $coverage_cutoff)
{
	print STDERR "warning: coverage-cutoff not set, assuming it is 1\n";
	$coverage_cutoff = 1;
}
elsif ($coverage_cutoff !~ /^\d+$/)
{
	die "coverage-cutoff=$coverage_cutoff is invalid\n".usage;
}

my %vcf_files;
my %mpileup_files;

my $dh;
# fill table vcf_files with entries like
#  vcf1 => dir/vcf1.vcf.gz
#  vcf2 => dir/vcf2.vcf.gz
opendir($dh, $vcf_dir) or die "error opening directory $vcf_dir: $!";
%vcf_files = map { /^(.*)\.vcf\.gz$/; $1 => "$vcf_dir/$_"} grep { /\.vcf\.gz$/ } readdir($dh);
closedir($dh);

die "No *.vcf.gz files found in $vcf_dir.  Perhas you need to compress and index with 'tabix' tools\n".
"Example: bgzip file.vcf; tabix -p vcf file.vcf.gz" if (keys(%vcf_files) <= 0);

my $total_samples = (keys %vcf_files);

# create table of mpileup files corresponding to input freebayes/variant vcf files
# assumes files are named with the same prefix
# ex vcf-dir/file1.vcf.gz and mpileup-dir/file1.vcf.gz
my $mpileup_table = create_mpileup_table(\%vcf_files, $vcf_dir, $mpileup_dir);
if (not defined $mpileup_table)
{
	die "Error: vcf-dir contains unmatched files in mpileup-dir";
}
else
{
	%mpileup_files = %{$mpileup_table};
}

# fill in variants for each vcf file
my $vcf_data = parse_variants(\%vcf_files,$requested_cpus);
my $mpileup_data = parse_mpileup(\%mpileup_files, $vcf_data,$requested_cpus);



my $invalid_pos;

$invalid_pos = parse_invalid($invalid) if $invalid;

my @samples_list = sort {$a cmp $b } keys %vcf_files;
my %chromosome_align;
my $unique_count = 1;
my %name_map; # used to map sample name to other information
my %sample_map; # keeps track of which samples have which unique ids (so we can properly increment unique_count)
my %total_positions_map; # keep track of total positions, and if valid/not
for my $chromosome (keys %$vcf_data)
{
	my ($alignment,$total_positions) = variants_alignment($vcf_data->{$chromosome}, $chromosome, $reference, \@samples_list, $mpileup_data, $coverage_cutoff,$invalid_pos,$keep_ambiguous);
	$total_positions_map{$chromosome} = $total_positions;
	for my $sample (sort {$a cmp $b} keys %$alignment)
	{
		next if (@{$alignment->{$sample}->{'positions'}} <= 0); # no alignments

		if (not exists $sample_map{$sample})
		{
			$sample_map{$sample} = "sample$unique_count";
			$unique_count++;
		}

		my $sample_id = $sample_map{$sample};
		my $sample_name = ($uniquify) ? $sample_id : $sample;

		if (not defined $chromosome_align{$sample})
		{
			$name_map{$sample_name} = "$sample\t$chromosome:".join('|',@{$alignment->{$sample}->{'positions'}});
			
			$chromosome_align{$sample} = {'header' => $sample_name};
			$chromosome_align{$sample}->{'data'} = $alignment->{$sample}->{'alignment'};
		}
		else
		{
			$name_map{$sample_name} .= " $chromosome:".join('|',@{$alignment->{$sample}->{'positions'}});

			$chromosome_align{$sample}->{'data'} .= $alignment->{$sample}->{'alignment'};
		}
	}
}
# print alignment
my $aln = Bio::SimpleAlign->new();
for my $sample (sort {$a cmp $b} keys %chromosome_align)
{
	my $id = $chromosome_align{$sample}->{'header'};
	my $data = $chromosome_align{$sample}->{'data'};

	if ($keep_ambiguous)
	{
		die "error: SNP alignment for $sample contains an invalid character"
			if ($data =~ /[^ATCGN]/);
	}
	else
	{
		die "error: SNP alignment for $sample contains an invalid character"
			if ($data =~ /[^ATCG]/);
	}

	my $seq = Bio::LocatableSeq->new(-seq => $data, -id => $id);
	$aln->add_seq($seq);
}

# sets displayname for each sequence
for my $seq_id ($aln->each_seq)
{
	my $start = $seq_id->start;
	my $end = $seq_id->end;
	my $id = $seq_id->id."/$start-$end";
	$aln->displayname($id, $seq_id->id);
}

# check if alignment is flush
die "Alignment blocks are not all of the same length" if (not $aln->is_flush());

for my $format (@formats)
{
	my $output_file = "$output_base.".$valid_formats{$format};
	my $io = Bio::AlignIO->new(-file => ">$output_file", -format => $format);
	$io->write_aln($aln);
	print STDERR "Alignment written to $output_file\n";
}

my $valid_positions = "$output_base-positions.tsv";
# print snp stats
print "# Command Line\n";
print "# $command_line\n";
print "# SNP statistics\n";
print "# Processed $total_samples samples\n";
print "# Total variant called SNPs processed: ".$snp_info->{'total'},"\n";
print "#\tRemoved ".$snp_info->{'removed'}->{'insertions'}," insertions\n";
print "#\tRemoved ".$snp_info->{'removed'}->{'deletions'}," deletions\n";
print "#\tRemoved ".$snp_info->{'removed'}->{'multi'}," multi\n";
print "#\tRemoved ".$snp_info->{'removed'}->{'other'}," other\n";
print "# Total Valid Positions: ".$snp_info->{'positions'},"\n";
print "# Total SNPs to process: ".($snp_info->{'positions'}*$total_samples)."\n";
print "#\tSNPs called as N's:\n";
print "#\t\tLow Coverage: ".$snp_info->{'snps'}->{'filtered-coverage'}."\n";
print "#\t\tVariant/mpileup differences: ".$snp_info->{'snps'}->{'filtered-mpileup'},"\n";
print "#\t\tInvalid Position based on user provided file: ".$snp_info->{'snps'}->{'filtered-invalid'},"\n" if $invalid_pos;
print "#\tValid SNPs for analysis: ".$snp_info->{'snps'}->{'kept'},"\n";
print "# Positions file in $valid_positions\n";

# print other information
print "#\n#AlnName\tSampleName\tPositions\n";
for my $name (sort keys %name_map)
{
	print "$name\t",$name_map{$name},"\n";
}

open(my $vfh, ">$valid_positions") or die "Could not open $valid_positions: $!";
print $vfh "#Chromosome\tPosition\tStatus\tReference\t";
my @samples_sorted_list = sort {$a cmp $b} @samples_list;
print $vfh join("\t",@samples_sorted_list);
print $vfh "\n";
for my $chr (keys %total_positions_map)
{
	my $positions = $total_positions_map{$chr};
	for my $pos (sort {$a <=> $b} keys %$positions)
	{
		my $samples = $positions->{$pos}->{'samples'};
		my $ref = $positions->{$pos}->{'ref'};
		print $vfh "$chr\t$pos\t".$positions->{$pos}->{'status'}."\t$ref\t";
		my $first = 1;
		die "error in total_positions_map, for $chr:$pos, not enough sample entries" if (@samples_sorted_list != scalar(keys %$samples));
		my $id = 0;
		for my $sample (sort {$a cmp $b } keys %$samples)
		{
			die "error: sample name $sample different from ".$samples_sorted_list[$id] if ($samples_sorted_list[$id] ne $sample);
			if ($first)
			{
				$first = 0;
				print $vfh $samples->{$sample};
			}
			else
			{
				print $vfh "\t".$samples->{$sample};
			}
			$id++;
		}
		print $vfh "\n";
	}
}
close($vfh);
