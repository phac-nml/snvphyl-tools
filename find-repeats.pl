#!/usr/bin/env perl
# filter-repeats.pl
# Purpose:  Given a nucmer repeats file, removes the "odd" repeats.

use warnings;
use strict;

use File::Basename;
use Cwd qw(abs_path getcwd);

my ($input) = @ARGV;
my $output = "output";
my $min_length;
my $min_pid;

sub usage
{
	"Usage: $0 [reference.fasta]\n".
	"Parameters:\n".
	"\t[reference.fasta]:  A fasta reference file to search for repeats.";
}

sub run_nucmer
{
	my ($input,$output_dir, $log_fh) = @_;

	my $input_base = basename($input, '.fasta');
	my $delta_prefix = "${input_base}_${input_base}";
	my $delta_file = "$output_dir/$delta_prefix.delta";

	my $command = "nucmer --maxmatch --nosimplify --prefix=$delta_prefix $input $input 1> $delta_prefix.out.log 2> $delta_prefix.err.log";

	my $cwd = getcwd;

	chdir $output_dir;
	print $log_fh "cd $output_dir\n";

	print $log_fh "$command\n";
	system($command) == 0 or die "Could not execute command '$command'";

	chdir $cwd;
	print $log_fh "cd $cwd\n\n";

	return $delta_file;
}

sub show_coords
{
	my ($delta_file,$min_length,$min_pid, $log_fh) = @_;

	my $coords_file = "$delta_file.coords";

	my $command = "show-coords -r -I $min_pid -L $min_length -TH $delta_file 1> $coords_file 2> $coords_file.err.log";

	print $log_fh "$command\n";

	system($command) == 0 or die "Could not execute command '$command'";

	return $coords_file;
}

# MAIN
$min_length = 150;
$min_pid = 90;
die "Error: no input file defined\n".usage if (not defined $input);
die "Error: file $input does not exist" if (not -e $input);
die "error: no output defined" if (not defined $output);
die "error: output directory $output already exists" if (-e $output);
mkdir $output if (not -e $output);

my $full_output = abs_path($output);
my $full_input = abs_path($input);
my $log = "$full_output/log.txt";
my $input_base = basename($input, '.fasta');

my $pos_good = "$full_output/$input_base.good.txt";
my $pos_bad = "$full_output/$input_base.bad.txt";

open(my $log_fh,">$log") or die "Could not open $log: $!";

my $delta_file = run_nucmer($full_input,$full_output,$log_fh);
my $coords_file = show_coords($delta_file,$min_length,$min_pid,$log_fh);

open(my $gfh, ">$pos_good") or die "Could not open $pos_good: $!";
open(my $bfh, ">$pos_bad") or die "Could not open $pos_bad: $!";

open(my $fh, "<$coords_file") or die "Could not open $coords_file: $!";
while(my $line = readline($fh))
{
	chomp $line;
	my @fields = split(/\t/,$line);
	die "Error: invalid show-coords file in $coords_file" if (@fields != 9);

	my ($start1,$end1,$start2,$end2,$length1,$length2,$pid,$contig1,$contig2) = @fields;
	my $diff_start = $start2-$start1;
	my $diff_end = $end2-$end1;

	if ($contig1 eq $contig2 and
		($diff_start == 0 or $diff_end == 0))
	{
		print $bfh "$line\n";
	}
	else
	{
		print $gfh "$line\n";
		print "$line\n";
	}
}
close($fh);
close($log_fh);
close($gfh);
close($bfh);
