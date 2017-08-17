#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use File::Temp 'tempfile';
use Getopt::Long;
use Template;
use File::Basename;

my $script_dir = $FindBin::Bin;
my $find_positions_used_bin = "$script_dir/../bin/find-positions-used.pl";
my $delete_temp = 1;
my $verbose = 1;

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
	if (defined readline($a_out_h))
	{
		$success = 0;
		fail("actual file $actual_out_file has more lines than expected file $expected_out_file");
	}
        close($a_out_h);

        return $success;
}

sub run_case
{
	my ($name,$reference_file,$core_gff_file,$bad_pos_file,$expected_out) = @_;

	my ($afh, $actual_out_file) = tempfile("find-used-pos.actual_out.XXXXXX", TMPDIR => 1, UNLINK => $delete_temp);
	close($afh);
	print STDERR "Actual Out: $actual_out_file\n" if ($verbose);

	my $command = "$find_positions_used_bin -c $core_gff_file -r $reference_file -b $bad_pos_file -t 1> $actual_out_file";
	$command .= " 2> /dev/null" if (not $verbose);
	print STDERR "$command\n" if ($verbose);
	system($command) == 0 or die "Could not execute command $command: $!";
	pass("pass $name") if (compare_files($expected_out, $actual_out_file));
}

### MAIN ###
my ($help);
my $usage = "Usage: $0 [-h|--help] [-v|--verbose]\n";
my $tests_dir = "$script_dir/find_positions_used_data";

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

opendir(my $test_in,$tests_dir) or die "Could not open $tests_dir: $!";
my @tests = map {"$tests_dir/$_"} sort {$a <=> $b} grep {$_ !~ /^\./} readdir($test_in);
closedir($test_in);

for my $dir (@tests)
{
	print "\nTesting $dir:\n" if ($verbose);

	my $reference_file = "$dir/input/reference.fasta";
	my $core_gff_file = "$dir/input/test.gff";
	my $bad_pos_file = "$dir/input/bad_pos.tsv";

	my $expected_out = "$dir/output/expected_out";

	run_case($dir,$reference_file,$core_gff_file,$bad_pos_file,$expected_out);
}

done_testing();
