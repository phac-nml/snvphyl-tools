#!/usr/bin/env perl
# check_snps_nucmer.pl
# Purpose:  Given a pseudoalign-positions.tsv and list of closed/finished fasta files
# and a reference, runs nucmer to verify the SNPs called from the reference mapping pipeline to those in pseudoalign-positions.tsv.

use warnings;
use strict;

use Getopt::Long;
use Set::Scalar;
use File::Basename;
use Bio::SeqIO;

use FindBin;
use lib $FindBin::Bin.'/lib';

use PositionsTable;
use CorePositions;
use NucmerPositionsChecker;

my ($input_align,$genome,$output,$reference,$bad_positions_file,$core_positions_file,$verbose);
$verbose = 0;
my $keep_temp = 1;

sub generate_genome_file_mapping
{
	my ($genomes_dir, $genome_names) = @_;

	my %genome_file_map; # maps {genome_name => closed_genome_file}
	opendir(my $d, $genomes_dir) or die "Could not open $genomes_dir: $!";
	while (my $genome_file = readdir $d)
	{
		next if ($genome_file !~ /\.fasta$/);
		my $curr_file = "$genomes_dir/$genome_file";
		print STDERR "curr_file=$curr_file\n";
		print STDERR "genome_file=$genome_file\n";

		# find genome name which appears at the beginning of the file name 
		for my $curr_name (@$genome_names)
		{
			print STDERR "checking $curr_name against $genome_file\n";
			if ($genome_file =~ /^$curr_name/)
			{
				if (defined $genome_file_map{$curr_name})
				{
					die "Error: found duplicate files for $curr_name: $curr_file and ".
						$genome_file_map{$curr_name}."\n";
				}
				else
				{
					$genome_file_map{$curr_name} = $curr_file;
				}
			}
		}
	}
	closedir($d);

	return \%genome_file_map;
}

sub generate_reference_contig_map
{
	my ($reference) = @_;

	my %contig_map;
	my $ref_io = Bio::SeqIO->new(-file=>"<$reference", -format=>"fasta");
	die "error: could not read $reference" if (not defined $ref_io);

	while (my $seq = $ref_io->next_seq)
	{
		$contig_map{$seq->display_id} = $seq;
	}

	return \%contig_map;
}

sub get_reference_length
{
	my ($reference) = @_;
	my $length = 0;

	my $io = Bio::SeqIO->new(-file=>"<$reference", -format=>"fasta");
	die "error: could not read $reference" if (not defined $io);

	while (my $seq = $io->next_seq)
	{
		$length += $seq->length;
	}

	return $length;
}

sub determine_genome_name
{
	my ($genomes_core_snp, $genome) = @_;
	my $genome_name = undef;
	my $genome_file = basename($genome, '.fasta');

	for my $curr_name (keys %$genomes_core_snp)
	{
		if ($genome_file =~ /^$curr_name/)
		{
			die "error: already exists name in table for genome $genome, name $genome_name" if (defined $genome_name);
			$genome_name = $curr_name;
		}
	}

	return $genome_name;
}

sub build_pipeline_set
{
	my ($genome_name, $genomes_core_snp) = @_;

	my $pipeline_set = Set::Scalar->new;

	my $genome_snps_set = $genomes_core_snp->{$genome_name};
	if (not defined $genome_snps_set)
	{
		warn "warning: no genome_snps_set for $genome_name defined";
	}

	for my $chrom (keys %$genome_snps_set)
	{
		my $pos_map = $genome_snps_set->{$chrom};
		for my $pos (keys %$pos_map)
		{
			my $ref = $pos_map->{$pos}->{'reference'};
			my $alt = $pos_map->{$pos}->{'alternative'};

			$pipeline_set->insert("$chrom\t$pos\t$ref\t$alt");
		}
	}

	return $pipeline_set;
}

sub print_snp_results
{
	my ($reference,$genome,$bad_positions_file,$core_positions_file,$pipeline_set,$nucmer_position_checker,$core_positions,$genome_core_snp_count) = @_;

	my $nucmer_set = $nucmer_position_checker->get_nucmer_snp_set;
	my $nucmer_ref_set = $nucmer_position_checker->get_nucmer_ref_set;
	my $nucmer_set_core_pos = $nucmer_position_checker->get_nucmer_snp_set_core_pos;
	my $nucmer_ref_set_core_pos = $nucmer_position_checker->get_nucmer_ref_set_core_pos;

	my $reference_base = basename($reference);
	my $genome_base = basename($genome);
	my $bad_positions_base = basename($bad_positions_file);
	my $core_positions_base = basename($core_positions_file);

	my $nucmer_all = ($nucmer_set + $nucmer_ref_set);
	my $nucmer_core_all = ($nucmer_set_core_pos + $nucmer_ref_set_core_pos);

	# compare sets of snps
	my $intersection = $pipeline_set * $nucmer_all;
	my $uniq_pipeline = $pipeline_set - $nucmer_all;
	my $uniq_nucmer = $nucmer_all - $pipeline_set;
	
	my $intersection_core_pos = $pipeline_set * $nucmer_core_all;
	my $uniq_pipeline_core = $pipeline_set - $nucmer_core_all;
	my $uniq_nucmer_core = $nucmer_core_all - $pipeline_set;
	
	my $total_bases_kept = scalar(keys %$core_positions);
	my $total_bases_reference = get_reference_length($reference);
	
	my $nucmer_snps = $nucmer_set->size;
	my $nucmer_filtered_snps = $nucmer_set_core_pos->size;
	
	print "Reference\tGenome\tCore Positions File\tBad Positions File\tTotal Reference Length\tTotal Length Used\t% Used\tCore Pipeline Positions\tCore Pipeline SNPs\tNucmer Positions\tNucmer SNPs\tNucmer Filtered Positions\tNucmer Filtered SNPs\tIntersection\tUnique Core Pipeline\tUnique Nucmer\t% True Positive\t% False Positive\t% False Negative\n";
	print "$reference_base\t$genome_base\t$core_positions_base\t$bad_positions_base\t$total_bases_reference\t$total_bases_kept\t";
	printf "%0.1f\t",($total_bases_kept/$total_bases_reference)*100;
	print $pipeline_set->size."\t$genome_core_snp_count\t".$nucmer_all->size."\t$nucmer_snps\t".$nucmer_core_all->size."\t$nucmer_filtered_snps\t".
		$intersection_core_pos->size."\t".$uniq_pipeline_core->size."\t".$uniq_nucmer_core->size."\t";
	
	my $true_positive;
	my $false_positive;
	my $false_negative;
	if ($nucmer_core_all->size > 0)
	{
		$false_negative = sprintf "%0.1f",($uniq_nucmer_core->size/$nucmer_core_all->size)*100;
	}
	else
	{
		$false_negative = 'undefined';
	}
	
	if ($pipeline_set->size > 0)
	{
		$false_positive = sprintf "%0.1f",($uniq_pipeline_core->size/$pipeline_set->size)*100;
		$true_positive = sprintf "%0.1f",($intersection_core_pos->size/$pipeline_set->size)*100;
	}
	else
	{
		$false_positive = 'undefined';
		$true_positive = 'undefined';
	}
	
	print "$true_positive\t$false_positive\t$false_negative\n";
	
	open(my $oh, ">$output") or die "Could not open file $output for writing: $!";
	print $oh "Working on $input_align\n";
	print $oh "Working with core positions $core_positions_file\n";
	print $oh "Working with bad positions $bad_positions_file\n";
	print $oh "Working with genome $genome\n";
	print $oh "Working with reference $reference\n\n";
	
	print $oh "Intersection\t".$intersection_core_pos->size."\n";
	print $oh "Contig\tPosition\t$reference_base\t$genome_base\n";
	for my $e (sort $intersection_core_pos->elements)
	{
		print $oh "$e\n";
	}
	
	print $oh "\nUnique to Nucmer\t".$uniq_nucmer_core->size."\n";
	print $oh "Contig\tPosition\t$reference_base\t$genome_base\n";
	for my $e (sort $uniq_nucmer_core->elements)
	{
		print $oh "$e\n";
	}
	
	print $oh "\nUnique to Core Pipeline\t".$uniq_pipeline_core->size."\n";
	print $oh "Contig\tPosition\t$reference_base\t$genome_base\n";
	for my $e (sort $uniq_pipeline->elements)
	{
		print $oh "$e\n";
	}
	
	close($oh);
}

sub usage
{
	"Usage: $0 -i [pseudoalign-positions.tsv] -r [reference file] -g [genome file] [-v]\n".
	"Parameters:\n".
	"\t-i|--input-align:  Input file (pseudoalign-positions.tsv generated by snp pipeline)\n".
	"\t-r|--reference: Reference genome.\n".
	"\t-o|--output: File to output detailed information\n".
	"\t-b|--bad-positions: Bad positions file\n".
	"\t-c|--core-positions: Core positions file\n".
	"\t-g|--genome:  Genome file to compare to.\n".
	"\t-v|--verbose\n";
}

# MAIN
if (!GetOptions('i|input-align=s' => \$input_align,
		'r|reference=s' => \$reference,
		'b|bad-positions=s' => \$bad_positions_file,
		'c|core-positions=s' => \$core_positions_file,
		'o|output=s' => \$output,
		'v|verbose' => \$verbose,
		'g|genome=s' => \$genome))
{
	die "Invalid option\n".usage;
}

die "Error: no input file defined\n".usage if (not defined $input_align);
die "Error: file $input_align does not exist\n".usage if (not -e $input_align);
die "Error: no genome file defined\n".usage if (not defined $genome);
die "Error: genome=$genome does not exist\n".usage if (not -e $genome);
die "Error: reference is not defined\n".usage if (not defined $reference);
die "Error: reference=$reference does not exist\n".usage if (not -e $reference);
die "Error: no bad-positions defind" if (not defined $bad_positions_file);
die "Error: bad-positions=$bad_positions_file does not exist" if (not -e $bad_positions_file);
die "Error: no core-positions defind" if (not defined $core_positions_file);
die "Error: core-positions=$core_positions_file does not exist" if (not -e $core_positions_file);
die "Error: no output file defined" if (not defined $output);


$keep_temp = 0 if ($verbose);

my $positions_table_parser = PositionsTable->new($verbose);
my ($genomes_core_snp,$genomes_core_snp_count) = $positions_table_parser->read_table($input_align);

my $genome_name = determine_genome_name($genomes_core_snp,$genome);
die "error: no entry in table $input_align for $genome" if (not defined $genome_name);

my $genome_core_snp_count = $genomes_core_snp_count->{$genome_name};
die "error: could not find SNP count in $input_align for $genome_name" if (not defined $genome_core_snp_count);

my $core_positions_parser = CorePositions->new($core_positions_file,$bad_positions_file,'tsv');
my $core_positions = $core_positions_parser->get_core_positions;

my $pipeline_set = build_pipeline_set($genome_name, $genomes_core_snp);

my $nucmer_position_checker = NucmerPositionsChecker->new($reference,$genome,$genome_name,$genomes_core_snp,$core_positions,$verbose);

print_snp_results($reference,$genome,$bad_positions_file,$core_positions_file,$pipeline_set,$nucmer_position_checker,$core_positions,$genome_core_snp_count);
