#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use File::Temp 'tempdir';
use File::Compare;

my $script_dir = $FindBin::Bin;

my $positions2snv_bin = "$script_dir/../positions2snv_alignment.pl";
my $positions2snv_dir = "$script_dir/positions2snv_alignment";

### MAIN

my $temp_dir = tempdir(CLEANUP => 0);
my $output_phylip = "$temp_dir/output.phy";
my $stdout_phylip = "$temp_dir/stdout-phylip";
my $input_positions = "$positions2snv_dir/empty_table/positions.tsv";

print "Testing case of empty SNV table\n";
my $command = "$positions2snv_bin -i $input_positions -o $output_phylip -f phylip > $stdout_phylip";
system($command) == 0 or die "Error executing $command\n";
ok(not (-e $output_phylip), "phylip file $output_phylip not produced for empty input table");
ok(`grep "No valid positions were found. Not creating empty alignment file" $stdout_phylip`, "Valid message on empty table");

opendir(my $in_h,$positions2snv_dir) or die "Could not open $positions2snv_dir: $!";
my @cases = sort {$a <=> $b} grep {$_ !~ /empty_table/} grep {$_ !~ /^\./} readdir($in_h);
closedir($in_h);

print "Testing all cases in $positions2snv_dir\n";
for my $case (@cases)
{
	my $case_dir = "$positions2snv_dir/$case";
	my $temp_dir = tempdir(CLEANUP => 0);

	print "\n### Testing $case_dir ###\n";

	my $done_testing = 0;
	my $input_positions = "$case_dir/positions.tsv";

	my $expected_phylip = "$case_dir/expected.phy";
	my $expected_fasta = "$case_dir/expected.fasta";

	my $output_phylip = "$temp_dir/output.phy";
	my $stdout_phylip = "$temp_dir/stdout-phylip";
	my $output_fasta = "$temp_dir/output.fasta";
	my $stdout_fasta = "$temp_dir/stdout-fasta";

	my $command = "$positions2snv_bin -i $input_positions -o $output_phylip -f phylip > $stdout_phylip";
	system($command) == 0 or die "Error executing $command\n";
	ok(compare($expected_phylip, $output_phylip) == 0, "phylip: $expected_phylip == $output_phylip");

	$command = "$positions2snv_bin -i $input_positions -o $output_fasta -f fasta > $stdout_fasta";
	system($command) == 0 or die "Error executing $command\n";
	ok(compare($expected_fasta, $output_fasta) == 0, "fasta: $expected_fasta == $output_fasta");
}

done_testing();
