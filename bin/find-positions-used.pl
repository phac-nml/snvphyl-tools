#!/usr/bin/env perl
# find-positions-used.pl
# Purpose:  Given a core.gff and bad_pos.tsv file, find the total
# number of positions used.

use warnings;
use strict;

use FindBin;
use Getopt::Long;
use Bio::SeqIO;

use lib $FindBin::Bin.'/lib';

use CorePositions;

my ($core_pos,$bad_pos,$reference,$verbose,$tab);
my %banned;
$verbose = 0;

sub get_reference_length
{
	my ($reference) = @_;
	my $length = 0;

	my $ref_i = Bio::SeqIO->new(-file=>"<$reference", -format=>"fasta");
	die "error: could not open $reference" if (not defined $ref_i);

	while (my $contig = $ref_i->next_seq)
	{
		$length += $contig->length;
	}

	return $length;
}

sub usage
{
	"Usage: $0 -c [core_pos.gff] -b [bad_pos.tsv] -r [reference.fasta] [-v]\n".
	"Parameters:\n".
	"\t-c|--core-pos:  Core positions file from vcf2core.pl (GFF format).\n".
	"\t-b|--bad-pos: Bad positions file, format of\n".
	"\t\tchromosome\tstart\tend\n".
	"\t-t|--tab: Print in tab-deliminted format\n".
	"\t-r|--reference:  Reference fasta file.\n";
}

# MAIN
if (!GetOptions('c|core-pos=s' => \$core_pos,
		'b|bad-pos=s' => \$bad_pos,
		't|tab' => \$tab,
		'r|reference=s' => \$reference))
{
	die "Invalid option\n".usage;
}

die "Error: no core-pos file defined\n".usage if (not defined $core_pos);
die "Error: file $core_pos does not exist" if (not -e $core_pos);
die "Error: no reference file defined\n".usage if (not defined $reference);
die "Error: file $reference does not exist" if (not -e $reference);
die "Error: no bad positions file defined\n".usage if (not defined $bad_pos);
die "Error: bad positions file $bad_pos does not exist\n".usage if (not -e $bad_pos);

my $core_positions = CorePositions->new($core_pos,$bad_pos,'gff');
my $ref_length = get_reference_length($reference);
my $core_count = $core_positions->get_core_pos_count;
my $used_count = $core_positions->get_used_pos_count;
my $percent_core = ($core_count/$ref_length)*100;
my $percent_used = ($used_count/$ref_length)*100;

if (not $tab)
{
	print "Date: ".`date`;
	print "Reference: $reference, length: $ref_length\n";
	print "Core Positions: $core_pos\n";
	print "Bad Positions: $bad_pos\n";
	printf "Core count: %i (%0.2f%%)\n",$core_count,$percent_core;
	printf "Used count: %i (%0.2f%%)\n",$used_count,$percent_used;
}
else
{
	print "Reference Length\tCore Length\tCore Percent\tUsed Length\tUsed Percent\n";
	printf "%i\t%i\t%0.2f\t%i\t%0.2f\n",$ref_length,$core_count,$percent_core,$used_count,$percent_used;
}
