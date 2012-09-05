#!/usr/bin/env perl
# vcf2pseudoalignmnet
# Purpose:  Given a set of *.vcf files, examines all called SNPs and generates a 'pseudoalignment'
#  of all the high quality SNPs/base calls (to be used for further phylogenetic analysis).

use warnings;
use strict;

use Getopt::Long;
use Storable 'dclone';

use Bio::AlignIO;
use Bio::SimpleAlign;
# VCF module from vcftools: http://vcftools.sourceforge.net/index.html
use Vcf;

my $verbose;
my $unknown_base = 'N';

my $snp_info = {'removed' => {'insertions' => 0, 'deletions' => 0, 'multi' => 0, 'other' => 0},
		'total' => 0,
		'positions' => 0,
		'snps' => {'filtered-coverage' => 0, 'filtered-mpileup' => 0, 'kept' => 0}
		};

sub usage
{
	"Usage: $0 --vcf-dir [vcf dir] --mpileup-dir [mpileup dir] --output [output alignment file]\n".
	"Parameters:\n".
	"\t--vcf-dir: The directory containing the vcf files.\n".
	"\t--mpileup-dir: Directory containing the vcf files produced by 'samtools mpileup \$file | bcftools view -cg' (*.vcf).\n".
	"\t-o|--output:  The output file for the alignment\n".
	"Options:\n".
	"\t-u|--uniquify:  Make the seq names unique and print mapping file to real names (for phylip format limitations)\n".
	"\t-r|--reference:  The name of the reference to use in the alignment (default: reference)\n".
	"\t-f|--format:  The format to output the alignment to, one of the Bio::AlignIO supported formats (default: fasta)\n".
	"\t-c|--coverage-cutoff:  The cutoff for coverage to include a reference base (default: 1)\n".
	"\t--verbose:  More information printed\n".
	"\t-h|--help:  Help\n";
}

sub check_subset_of
{
	my ($set1, $set1_dir, $set2, $set2_dir) = @_;

	for my $key1 (keys %$set1)
	{
		if (not defined $set2->{$key1})
		{
			die "Could not find file in $set2_dir matching file ".$set1->{$key1}." in $set1_dir";
		}
		print STDERR $set1->{$key1}." matches ".$set2->{$key1}."\n" if ($verbose);
	}
}

sub variant_info_to_hash
{
	my ($info_string) = @_;

	my @info = split(/;/,$info_string);
	my %info_hash = map {my @a=split(/=/,$_); $a[0] => ((@a >= 1) ? $a[1] : undef)} @info;

	return \%info_hash;
}

# Parses position data structure (for a single chromosome) and builds alignment.
# Input:  positions_hash  The hash of all positions of variants
#	  samples_list  A list of all sample names to proccess
# Output:  An alignment of all valid variant positions (in FASTA format)
sub variants_alignment
{
	my ($positions_hash, $chromosome, $reference, $samples_list, $mpileup_data, $coverage_cutoff) = @_;
	my $alignment_string = undef;

	# stores pseudoalignment in form of
	# {sample => {alignment => 'string', positions => [pos array]}}
	my %alignment;
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
				}
				else
				{
					my $ref = $pileup_vcf->[3];
					my $alt = $pileup_vcf->[4];
					my %info_hash = %{variant_info_to_hash($pileup_vcf->[7])};

					my $coverage = $info_hash{'DP'};
					if (not defined $coverage)
					{
						die "Error: coverage not defined for sample $sample, vcf-line: ".join("\t",@$pileup_vcf);
					}
					elsif ($coverage < $coverage_cutoff)
					{
						print STDERR "fail for $sample:$chromosome:$pos, coverage=$coverage < cutoff=$coverage_cutoff\n" if ($verbose);
						$alignment_local->{$sample} = {'base' => $unknown_base, 'position' => $pos};
						$snp_info->{'snps'}->{'filtered-coverage'}++;
					}
					elsif ($alt ne '.')
					{
						print STDERR "fail for $sample:$chromosome:$pos, mpileup data gives alt=$alt (ref=$ref), but no variant called from vcf data\n" if ($verbose);
						$alignment_local->{$sample} = {'base' => $unknown_base, 'position' => $pos};
						$snp_info->{'snps'}->{'filtered-mpileup'}++;
					}
					else
					{
						print STDERR "pass for $sample:$chromosome:$pos, coverage=$coverage > cutoff=$coverage_cutoff and no variant called from mpileup data\n" if ($verbose);
						$alignment_local->{$sample} = {'base' => $ref_base, 'position' => $pos};

						$snp_info->{'snps'}->{'kept'}++;
					}
				}
			}
			else
			{
				$alignment_local->{$sample} = {'base' => $sample_hash->{$sample}->{'alt'}, 'position' => $pos};
				$snp_info->{'snps'}->{'kept'}++;
			}
		}

		# fill in overall data structure
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


	return \%alignment;
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
	my ($vcf_files) = @_;

	my %vcf_data;
	for my $vcf_key (keys %$vcf_files)
	{
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
					print STDERR "multi variant region found in $vcf_file, skipping: \"".join(' ',@$data),"\"\n" if ($verbose);
					$snp_info->{'removed'}->{'multi'}++;
				}
				elsif (length($ref) > 1 and length($ref) < length($alt))
				{
					print STDERR "length($ref) > 1 and length($ref) < length($alt) in $vcf_file, skipping: \"".join(' ',@$data),"\"\n" if ($verbose);
					$snp_info->{'removed'}->{'other'}++;
				}
				elsif (length($ref) > 1) # deletion in query strain
				{
					$snp_info->{'removed'}->{'deletions'}++;
					print STDERR "deletion found in $vcf_file, skipping: \"".join(' ',@$data),"\"\n" if ($verbose);
				}
				elsif (length($alt) > 1) # insertion in query strain
				{
					$snp_info->{'removed'}->{'insertions'}++;
					print STDERR "insertion found in $vcf_file, skipping: \"".join(' ',@$data),"\"\n" if ($verbose);
				}
				elsif((length($alt) == 1) and length($ref) == 1) # SNP
				{
					print STDERR "for file $vcf_file, keeping variant \"".join(' ',@$data),"\"\n" if ($verbose);
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
					print STDERR "invalid lengths for ref and query found in $vcf_file, skipping: \"".join(' ',@$data),"\"\n" if ($verbose);
					$snp_info->{'removed'}->{'other'}++;
				}

			}
		}
	}

	return \%vcf_data;
}

sub parse_mpileup
{
	my ($mpileup_files, $vcf_data) = @_;
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

	# read through each vcf file
	# place into mpileup_hash, keeping only vcf lines necessary (in sorted_positions)
	for my $sample (keys %$mpileup_files)
	{
		my $working_mpileup = dclone($positions_hash);
		$sample_pileup_hash{$sample} = $working_mpileup;

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
					if (exists $info_hash{'INDEL'})
					{
						print STDERR "skipping INDEL vcf-line found for $sample: '".join(" ",@$data)."\n";
					}
					elsif (defined $working_mpileup->{$chrom}->{$position})
					{
						die "Duplicate position found in $sample.\nvcf-line: ".join(" ",@$data);
					}
					else
					{
						print STDERR "mpileup info found for $sample:$chrom:$position\n" if ($verbose);
						$working_mpileup->{$chrom}->{$position} = $data;
					}
				}
			}
		}
	}

	return \%sample_pileup_hash;
}

############
### MAIN ###
############

my %valid_formats = map {$_ => 1} ('fasta', 'phylip', 'clustalw');

my $vcf_dir;
my $mpileup_dir;
my $output;
my $format;
my $reference;
my $coverage_cutoff;
my $uniquify;
my $help;

my $command_line = join(' ',@ARGV);

if (!GetOptions('vcf-dir|d=s' => \$vcf_dir,
		'mpileup-dir|b=s' => \$mpileup_dir,
		'format|f=s' => \$format,
		'output|o=s' => \$output,
		'reference|r=s' => \$reference,
		'coverage-cutoff|c=i' => \$coverage_cutoff,
		'uniquify|u' => \$uniquify,
		'help|h' => \$help,
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

die "output file undefined\n".usage if (not defined $output);

if (not defined $reference)
{
	print STDERR "reference name not defined, calling it 'reference'\n";
	$reference = 'reference';
}

$uniquify = 0 if (not defined $uniquify);

if (not defined $format)
{
	print STDERR "warning: format not defined, assuming fasta\n";
	$format = "fasta";
}
elsif (not defined $valid_formats{$format})
{
	die "unrecognized format '$format', must be one of '".join(' ', keys %valid_formats),"'\n";
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

# fill table depth_files with entries like
#  sample1 => dir/vcf1.vcf.gz
#  sample2 => dir/vcf2.vcf.gz
opendir($dh, $mpileup_dir) or die "error opening directory $mpileup_dir: $!";
%mpileup_files = map { /^(.*)\.vcf\.gz$/; $1 => "$mpileup_dir/$_"} grep { /\.vcf\.gz$/ } readdir($dh);
closedir($dh);

# check to make sure every vcf file has a corresponding depth file
# assumes files are named with the same prefix and just have the .depth or .vcf suffix changed
# ex file1.vcf and file1.depth
check_subset_of(\%vcf_files, $vcf_dir, \%mpileup_files, $mpileup_dir);

# fill in variants for each vcf file
my $vcf_data = parse_variants(\%vcf_files);
my $mpileup_data = parse_mpileup(\%mpileup_files, $vcf_data);

my @samples_list = keys %vcf_files;
my %chromosome_align;
my $unique_count = 1;
my %name_map; # used to map sample name to other information
my %sample_map; # keeps track of which samples have which unique ids (so we can properly increment unique_count)
for my $chromosome (keys %$vcf_data)
{
	my $alignment = variants_alignment($vcf_data->{$chromosome}, $chromosome, $reference, \@samples_list, $mpileup_data, $coverage_cutoff);
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

my $io = Bio::AlignIO->new(-file => ">$output", -format => $format);
$io->write_aln($aln);
print STDERR "Alignment written to $output\n";

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
print "#\tValid SNPs for analysis: ".$snp_info->{'snps'}->{'kept'},"\n";

# print other information
print "#\n#AlnName\tSampleName\tPositions\n";
for my $name (sort keys %name_map)
{
	print "$name\t",$name_map{$name},"\n";
}

