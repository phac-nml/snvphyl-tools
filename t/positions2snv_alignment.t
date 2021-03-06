#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use File::Temp 'tempdir';
use Bio::AlignIO;

my $script_dir = $FindBin::Bin;

my $positions2snv_bin = "$script_dir/../bin/positions2snv_alignment.pl";
my $positions2snv_dir = "$script_dir/snv_alignment_data";

sub compare_alignments
{
	my ($a, $b) = @_;

	my $aio = Bio::AlignIO->new(-file=>"<$a",-longid=>1);
	my $bio = Bio::AlignIO->new(-file=>"<$b",-longid=>1);

	my $a_aln = $aio->next_aln;
	my $b_aln = $bio->next_aln;

	#return 0 if ((defined $aio->next_aln) or (defined $bio->next_aln));
	return 0 if ($a_aln->num_sequences != $b_aln->num_sequences);
	
	for (my $i = 1; $i <= $a_aln->num_sequences; $i++)
	{
		my $aseq = $a_aln->get_seq_by_pos($i);
		my $bseq = $b_aln->get_seq_by_pos($i);

		return 0 if ($aseq->display_id ne $bseq->display_id);
		return 0 if ($aseq->seq ne $bseq->seq);
	}

	return 1;
}

### MAIN

my $temp_dir = tempdir(CLEANUP => 1);
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
	my $temp_dir = tempdir(CLEANUP => 1);

	print "\n### Testing $case_dir ###\n";

	my $done_testing = 0;
	my $input_positions = "$case_dir/positions.tsv";
	my $input_reference = "$case_dir/reference.fasta";

	my $expected_phylip = "$case_dir/expected.phy";
	my $expected_fasta = "$case_dir/expected.fasta";
	my $expected_fasta_all = "$case_dir/expected-all.fasta";
	my $expected_stdout = "$case_dir/expected-stdout-short";
	my $expected_stdout_all = "$case_dir/expected-stdout-short-all";

	my $output_phylip = "$temp_dir/output.phy";
	my $stdout_phylip = "$temp_dir/stdout-phylip";
	my $output_fasta = "$temp_dir/output.fasta";
	my $output_fasta_all = "$temp_dir/output-all.fasta";
	my $stdout_fasta = "$temp_dir/stdout-fasta";
	my $stdout_fasta_all = "$temp_dir/stdout-fasta-all";

	my $command = "$positions2snv_bin -i $input_positions -o $output_phylip -f phylip > $stdout_phylip";
	system($command) == 0 or die "Error executing $command\n";
	ok(compare_alignments($expected_phylip, $output_phylip), "phylip: $expected_phylip == $output_phylip");
	ok(`grep -f $expected_stdout $stdout_phylip -c | tr -d '\n'` == `wc -l $expected_stdout | cut -d ' ' -f 1 | tr -d '\n'`, "Every line in $expected_stdout matched actual output $stdout_phylip");

	$command = "$positions2snv_bin -i $input_positions -o $output_fasta -f fasta > $stdout_fasta";
	system($command) == 0 or die "Error executing $command\n";
	ok(compare_alignments($expected_fasta, $output_fasta), "fasta: $expected_fasta == $output_fasta");
	ok(`grep -f $expected_stdout $stdout_fasta -c | tr -d '\n'` == `wc -l $expected_stdout | cut -d ' ' -f 1 | tr -d '\n'`, "Every line in $expected_stdout matched actual output $stdout_fasta");

	$command = "$positions2snv_bin -i $input_positions --keep-all -o $output_fasta_all -f fasta > $stdout_fasta_all";
	system($command) == 0 or die "Error executing $command\n";
	ok(compare_alignments($expected_fasta_all, $output_fasta_all), "fasta: $expected_fasta_all == $output_fasta_all");
	ok(`grep -f $expected_stdout_all $stdout_fasta_all -c | tr -d '\n'` == `wc -l $expected_stdout_all | cut -d ' ' -f 1 | tr -d '\n'`, "Every line in $expected_stdout_all matched actual output $stdout_fasta_all");
}

done_testing();
