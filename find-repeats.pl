#!/usr/bin/env perl
# find-repeats.pl
# Purpose:  Runs numcer to find repeats and does some basic filtering.

use warnings;
use strict;

use File::Basename;
use Getopt::Long;
use File::Temp 'tempdir';
use Cwd qw(abs_path getcwd);

my $min_length_default = 150;
my $min_pid_default = 90;
my $min_length;
my $min_pid;
my $keep_temp = 0;

sub usage
{
	"Usage: $0 [reference.fasta] --min-length [length] --min-pid [pid]\n".
	"Parameters:\n".
	"\t[reference.fasta]:  A fasta reference file to search for repeats.\n".
	"Options:\n".
	"\t-l|--min-length: Minimum length of repeat region ($min_length_default).\n".
	"\t-p|--min-pid: Minimum PID of repeat region ($min_pid_default).\n".
	"\t-k|--keep-temp: Keep around temporary nucmer/coords files (no).\n";
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
if (!GetOptions('l|min-length=i' => \$min_length,
		'p|min-pid=i' => \$min_pid,
		'k|keep-temp' => \$keep_temp))
{
	die "Invalid option\n".usage;
}

my ($input) = @ARGV;
die "error: no input file defined\n".usage if (not defined $input);
die "error: file $input does not exist\n".usage if (not -e $input);

$min_length = $min_length_default if (not defined $min_length);
$min_pid = $min_pid_default if (not defined $min_pid);

die "error: $min_length must be positive\n".usage if ($min_length <= 0);
die "error: $min_pid must be between 0-100\n".usage if ($min_pid < 0 or $min_pid > 100);

print STDERR "Finding repeats on $input\n";
print STDERR "Using min-length $min_length, min-pid $min_pid\n";

my $output = tempdir('find_repeats_XXXXXX', TMPDIR => 0, CLEANUP => (not $keep_temp));

my $full_output = abs_path($output);
my $full_input = abs_path($input);
my $log = "$full_output/log.txt";
my $input_base = basename($input, '.fasta');

my $pos_good = "$full_output/$input_base.good.txt";
my $pos_bad = "$full_output/$input_base.bad.txt";

open(my $log_fh,">$log") or die "Could not open $log: $!";

print STDERR "Temporary directory: $full_output\n" if ($keep_temp);

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

		print "$contig1\t$start1\t$end1\n";
		print "$contig2\t$start2\t$end2\n";
	}
}
close($fh);
close($log_fh);
close($gfh);
close($bfh);
