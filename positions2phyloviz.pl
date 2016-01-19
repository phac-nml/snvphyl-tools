#!/usr/bin/env perl
# positions2phyloviz
# Purpose:  Given a snv_align-positions.tsv, converts to a file that can be imported into phyloviz

use warnings;
use strict;

use Getopt::Long;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;

my ($input,$output,$reference_name,$verbose);
$verbose = 0;

sub usage
{
	"Usage: $0 -i [snv_align-positions.tsv] -o [snv_align out] -f [alignment format] [-v]\n".
	"Parameters:\n".
	"\t-i|--input:  Input file (snv_align-positions.tsv generated by snv pipeline)\n".
	"\t-b|--output-base:  Output file base name\n".
	"\t--reference-name: Use passed name instead of default for reference\n".
	"\t--verbose: Print more information\n";

}

if (!GetOptions('i|input=s' => \$input,
		'b|output-base=s' => \$output,
		'reference-name=s' => \$reference_name,
		'v|verbose' => \$verbose))
{
	die "Invalid option\n".usage;
}

die "Error: no input file defined\n".usage if (not defined $input);
die "Error: file $input does not exist" if (not -e $input);
die "Error: no output file defined\n".usage if (not defined $output);

print "Date: ".`date`;
print "Working on $input\n";

open(my $fh, "<$input") or die "Could not open $input: $!";

my $line = readline($fh);
chomp($line);

die "Error: no header line defined in $input" if ($line !~ /^#Chromosome\tPosition\tStatus\tReference/);
my @values = split(/\t/,$line);
my (undef,undef,undef,@strains) = @values;
die "Error: no strains defined in $input" if (@strains <= 0);

# replace reference name
if ($strains[0] eq 'Reference' and defined $reference_name)
{
	$strains[0] = $reference_name;
}

# below arrays are parallel with @strains

# stores array of nucleotide data for each strain, parallel with @strains
my @nucleotide_data;

# array of positions processed, parallel to each element in @nucleotide_data
my @position_information;

my $valid_count=0;
my $invalid_count=0;
while($line = readline($fh))
{
	chomp $line;
	@values = split(/\t/,$line);

	my ($chrom,$pos,$status,@dna) = @values;

	if (scalar(@dna) != scalar(@strains))
	{
		die "Error: line $line does not have same number of entries as header for $input";
	}
	elsif ($status ne 'valid')
	{
		print STDERR "Skipping invalid line: $line\n" if ($verbose);
		$invalid_count++;
	}
	else
	{
		$valid_count++;

		# fill in position information
		push(@position_information, "${chrom}_${pos}");

		for (my $i = 0; $i < @dna; $i++)
		{
			die "error: DNA base (".$dna[$i].") not valid, must be [ACTG] for line \"$line\"" if ($dna[$i] !~ /[ACTG]/i);
			if (not defined $nucleotide_data[$i])
			{
				$nucleotide_data[$i] = $dna[$i];
			}
			else
			{
				$nucleotide_data[$i] .= "\t".$dna[$i];
			}
		}
	}
}
close($fh);

my $out_profile = "$output-profile.tsv";
my $out_strains = "$output-strains.tsv";

open(my $profile_h, ">$out_profile") or die "Could not open $out_profile for writing: $!";
open(my $strain_h, ">$out_strains") or die "Could not open $out_strains for writing: $!";
print $strain_h "ST\tStrain\n";
print $profile_h "ST\t".join("\t",@position_information)."\n";

# generate and print out sequence types, check which types (set of SNVs for a strain) are identical
# sequence types stored in table as {SNV profile => [initial_strain, type number]}
my %types;
my $curr_unused_type = 1; # start snv profile types at 1
my $duplicate_profiles = 0;
for (my $curr = 0; $curr < @strains; $curr++)
{
	my $strain_id = $strains[$curr];
	my $snv_profile = $nucleotide_data[$curr];
	my $curr_type;

	if (not exists $types{$snv_profile})
	{
		$curr_type = $curr_unused_type;
		$curr_unused_type++;

		$types{$snv_profile} = [$strain_id, $curr_type];

		# only print profile information if this is a unique sequence type
		print $profile_h "$curr_type\t$snv_profile\n";
	}
	else
	{
		my $previous_strain_id;
		($previous_strain_id,$curr_type) = @{$types{$snv_profile}};
		print STDERR "profile for $strain_id duplicate of $previous_strain_id, assigning to $curr_type\n" if ($verbose);
		$duplicate_profiles++;
	}

	# print information for strain
	print $strain_h "$curr_type\t$strain_id\n";
}
close($profile_h);
close($strain_h);

print "Kept $valid_count valid positions\n";
print "Skipped $invalid_count positions\n";
print "Found $duplicate_profiles duplicate SNV profiles\n";
print "SNV profile information written to $out_profile\n";
print "SNV type/strain mapping written to $out_strains\n";
