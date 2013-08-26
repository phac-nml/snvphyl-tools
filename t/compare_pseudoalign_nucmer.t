#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use File::Temp 'tempfile';
use Getopt::Long;

my $script_dir = $FindBin::Bin;
my $compare_snps_bin = "$script_dir/../compare_pseudoalign_nucmer.pl";
my $delete_temp = 1;
my $verbose = 0;

sub compare_files
{
        my ($expected_out_file,$actual_out_file) = @_;

        my $success = 1;

        open(my $out_h, $expected_out_file) or die "Could not open $expected_out_file: $!";
        open(my $a_out_h, $actual_out_file) or die "Could not open $actual_out_file: $!";
        while($success and (defined (my $expected_line = readline($out_h))))
        {
                my $actual_line = readline($a_out_h);
                if (not defined $actual_line)
                {
                        $success = 0;
                        fail("expected file $expected_out_file has more lines than actual file $actual_out_file");
                        next;
                }
                else
                {
                        chomp $expected_line;
                        chomp $actual_line;
                        if ($actual_line ne $expected_line)
                        {
                                is($actual_line,$expected_line,"lines \"$actual_line\" and \"$expected_line\" differ");
                                $success = 0;
                        }
                }
        }
        close($out_h);
        close($a_out_h);

        return $success;
}

sub build_input_files
{
	my ($reference, $fasta, $fasta_file_prefix, $pseudoalign) = @_;

	my ($rfh,$reference_file) = tempfile('compare_nucmer.reference.test.XXXXXX', TMPDIR => 1, UNLINK => $delete_temp);
	print $rfh $reference;
	close($rfh);
	print STDERR "Reference: $reference_file\n" if ($verbose);

	my ($ffh, $fasta_file) = tempfile("$fasta_file_prefix.compare_nucmer.XXXXXX", TMPDIR => 1, UNLINK => $delete_temp);
	print $ffh $fasta;
	close($ffh);
	print STDERR "FASTA: $fasta_file\n" if ($verbose);

	my ($pfh, $pseudoalign_file) = tempfile("compare_nucmer.pseudoalign.XXXXXX", TMPDIR => 1, UNLINK => $delete_temp);
	print $pfh $pseudoalign;
	close($pfh);
	print STDERR "Pseudoalign: $pseudoalign_file\n" if ($verbose);

	return ($reference_file, $fasta_file, $pseudoalign_file);
}

sub run_case
{
	my ($name, $reference, $fasta, $pseudoalign, $expected_out) = @_;

	my ($afh, $actual_out_file) = tempfile("compare_nucmer.actual_out.XXXXXX", TMPDIR => 1, UNLINK => $delete_temp);
	close($afh);
	print STDERR "Actual Out: $actual_out_file\n" if ($verbose);
	
	my $command = "$compare_snps_bin -i $pseudoalign -r $reference -g $fasta 1> $actual_out_file";
	print STDERR "$command\n" if ($verbose);
	system($command) == 0 or die "Could not execute command $command: $!";
	pass("pass $name") if (compare_files($expected_out, $actual_out_file));
}

sub build_expected_out_1
{
	my ($reference_file, $fasta_file) = @_;
	my $expected_out_1 =
	"#Reference\tGenome\tPipeline\tNucmer\tIntersection\tUniqPipeline\tUniqNucmer\tTruePositive\tFalsePositive\tFalseNegative\n".
	"$reference_file\t$fasta_file\t1\t1\t1\t0\t0\t1.000\t0.000\t0.000\n";

	my ($efh, $expected_out_file) = tempfile("compare_nucmer.expected_out.XXXXXX", TMPDIR => 1, UNLINK => $delete_temp);
	print $efh $expected_out_1;
	close($efh);
	print STDERR "Expected Out 1: $expected_out_file\n" if ($verbose);

	return $expected_out_file;
}

my $reference =
">chr\n".
"TGAAATCGAATCGGATTCG\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"TTTTTTTTTTTTTTT\n";

my $fasta_in_A =
">chr\n".
"TGAATTCGAATCGGATTCG\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"AAAAAAAAAAAAAAA\n".
"TTTTTTTTTTTTTTT\n";

my $fasta_in_B = 
">chr\n".
"ATCGATCGATCG\n";

my $example_pseudoalign =
"#Chromosome\tPosition\tStatus\tReference\tA\tB\n".
"chr\t5\tvalid\tA\tT\tA\n";

### MAIN ###
my ($help);
my $usage = "Usage: $0 [-h|--help] [-v|--verbose]\n";
my ($reference_file, $fasta_file, $pseudoalign_file, $expected_out_file);

if (!GetOptions('v|verbose' => \$verbose,
                'h|help' => \$help))
{
        die "Invalid option\n$usage";
}

if ($help)
{
	print $usage;
	exit 0;
}

if (not defined $verbose)
{
	$verbose = 0;
}
elsif ($verbose)
{
	$delete_temp = 0;
}

print "Testing $compare_snps_bin\n";

($reference_file, $fasta_file, $pseudoalign_file) = build_input_files($reference, $fasta_in_A, "A", $example_pseudoalign);
$expected_out_file = build_expected_out_1($reference_file, $fasta_file);
run_case("Case 1", $reference_file, $fasta_file, $pseudoalign_file, $expected_out_file);

done_testing();
