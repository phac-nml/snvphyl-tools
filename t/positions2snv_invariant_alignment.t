#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use File::Temp 'tempdir';
use File::Compare;
use Bio::SeqIO;

my $script_dir = $FindBin::Bin;

my $positions2snv_invariant_bin = "$script_dir/../positions2snv_invariant_alignment.pl";
my $positions2snv_dir = "$script_dir/snv_alignment_data";

### MAIN

my $temp_dir = tempdir(CLEANUP => 1);
my $output_invariant_dir = "$temp_dir/output";
my $output_phylip = "$output_invariant_dir/ref.phylip";
my $stdout_phylip = "$temp_dir/stdout-phylip";
my $input_positions = "$positions2snv_dir/empty_table/positions.tsv";
my $reference_file = "$positions2snv_dir/empty_table/reference.fasta";

print "Testing case of empty SNV table\n";
my $command = "$positions2snv_invariant_bin -i $input_positions -o $output_invariant_dir --reference-file $reference_file -f phylip > $stdout_phylip";
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
	my $temp_dir = tempdir(CLEANUP => 1);

	print "\n### Testing $case_dir ###\n";

	my $done_testing = 0;
	my $input_positions = "$case_dir/positions.tsv";
	my $input_reference = "$case_dir/reference.fasta";
	my $seq_io = Bio::SeqIO->new(-file=>"<$input_reference", -format=>"fasta");
	my @reference_ids;
	while (my $seq = $seq_io->next_seq)
	{
		push(@reference_ids, $seq->display_id.".fasta");
	}
	die "Error, no valid reference ids in file $input_reference for test" if (scalar(@reference_ids) == 0);

	my $expected_fasta_invariant_dir = "$case_dir/invariant";
	my $expected_fasta_invariant_merged_dir = "$case_dir/invariant-merged";
	my $expected_fasta_invariant_all_dir = "$case_dir/invariant-all";
	my $expected_fasta_invariant_merged_all_dir = "$case_dir/invariant-merged-all";
	my $expected_stdout_invariant = "$case_dir/expected-stdout-short";
	my $expected_stdout_invariant_merged = "$case_dir/expected-stdout-short-merged";
	my $expected_stdout_invariant_all = "$case_dir/expected-stdout-short-all";

	my $output_fasta_invariant_dir = "$temp_dir/invariant";
	my $output_fasta_invariant_merged_dir = "$temp_dir/invariant-merged";
	my $output_fasta_invariant_all_dir = "$temp_dir/invariant-all";
	my $output_fasta_invariant_merged_all_dir = "$temp_dir/invariant-merged-all";
	my $stdout_fasta_invariant = "$temp_dir/stdout-fasta-invariant";
	my $stdout_fasta_invariant_merged = "$temp_dir/stdout-fasta-invariant-merged";
	my $stdout_fasta_invariant_all = "$temp_dir/stdout-fasta-invariant-all";
	my $stdout_fasta_invariant_merged_all = "$temp_dir/stdout-fasta-invariant-merged-all";

	# Test writing alignment with only valid positions
	$command = "$positions2snv_invariant_bin -i $input_positions -o $output_fasta_invariant_dir -f fasta --reference-file $input_reference > $stdout_fasta_invariant";
	system($command) == 0 or die "Error executing $command\n";
	for my $ref_id (@reference_ids)
	{
		ok(compare("$expected_fasta_invariant_dir/$ref_id", "$output_fasta_invariant_dir/$ref_id") == 0, "fasta: $expected_fasta_invariant_dir/$ref_id == $output_fasta_invariant_dir/$ref_id");
	}
	ok(`grep -f $expected_stdout_invariant $stdout_fasta_invariant -c | tr -d '\n'` == `wc -l $expected_stdout_invariant | cut -d ' ' -f 1 | tr -d '\n'`, "Every line in $expected_stdout_invariant matched actual output $stdout_fasta_invariant");

	# Test alignment with both valid and invald positions
	$command = "$positions2snv_invariant_bin --keep-all -i $input_positions -o $output_fasta_invariant_all_dir -f fasta --reference-file $input_reference > $stdout_fasta_invariant_all";
	system($command) == 0 or die "Error executing $command\n";
	for my $ref_id (@reference_ids)
	{
		ok(compare("$expected_fasta_invariant_all_dir/$ref_id", "$output_fasta_invariant_all_dir/$ref_id") == 0, "fasta: $expected_fasta_invariant_all_dir/$ref_id == $output_fasta_invariant_all_dir/$ref_id");
	}
	ok(`grep -f $expected_stdout_invariant_all $stdout_fasta_invariant_all -c | tr -d '\n'` == `wc -l $expected_stdout_invariant_all | cut -d ' ' -f 1 | tr -d '\n'`, "Every line in $expected_stdout_invariant_all matched actual output $stdout_fasta_invariant_all");

	# Test merging alignment together into a single file (valid positions)
	$command = "$positions2snv_invariant_bin --merge-alignment -i $input_positions -o $output_fasta_invariant_merged_dir -f fasta --reference-file $input_reference > $stdout_fasta_invariant_merged";
	system($command) == 0 or die "Error executing $command\n";
	ok(compare("$expected_fasta_invariant_merged_dir/alignment_merged.fasta", "$output_fasta_invariant_merged_dir/alignment_merged.fasta") == 0, "fasta: $expected_fasta_invariant_merged_dir/alignment_merged.fasta == $output_fasta_invariant_merged_dir/alignment_merged.fasta");

	# Test merging alignment together into a single file (all positions)
	$command = "$positions2snv_invariant_bin --keep-all --merge-alignment -i $input_positions -o $output_fasta_invariant_merged_all_dir -f fasta --reference-file $input_reference > $stdout_fasta_invariant_merged_all";
	system($command) == 0 or die "Error executing $command\n";
	ok(compare("$expected_fasta_invariant_merged_all_dir/alignment_merged.fasta", "$output_fasta_invariant_merged_all_dir/alignment_merged.fasta") == 0, "fasta: $expected_fasta_invariant_merged_all_dir/alignment_merged.fasta == $output_fasta_invariant_merged_all_dir/alignment_merged.fasta");
}

done_testing();
