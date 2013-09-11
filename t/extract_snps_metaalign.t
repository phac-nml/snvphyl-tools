#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use File::Temp 'tempfile';
use Getopt::Long;

my $script_dir = $FindBin::Bin;
my $extract_snps_bin = "$script_dir/../extract_snps_metaalign.pl";
my $delete_temp = 1;

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

sub run_case
{
	my ($name, $example_file, $case_text, $command_line) = @_;

	my ($cfh,$case_expected) = tempfile('extract_snps_metaalign.test.XXXXXX', TMPDIR => 1, UNLINK => $delete_temp);
	print $cfh $case_text;
	close($cfh);
	
	my ($cfha, $case_actual) = tempfile('extract_snps_metaalign.test.XXXXXX', TMPDIR => 1, UNLINK => $delete_temp);
	close($cfha);
	
	my $command = "$extract_snps_bin -i $example_file -o $case_actual $command_line > /dev/null 2> /dev/null";
	#print "$command\n";
	system($command) == 0 or die "Could not execute command $command: $!";
	pass("pass $name") if (compare_files($case_expected, $case_actual));
}

my $example =
"#Chromosome\tPosition\tStatus\tReference\tA\tB\tC\n".
"chr\t1\tvalid\tA\tA\tA\tT\n".
"chr\t2\tvalid\tA\tC\tC\tT\n".
"chr\t3\tvalid\tA\tC\tG\tT\n".
"chr\t4\tfiltered-mpileup\tC\tC\tN\tT\n";

my $case1 = 
"#Chromosome\tPosition\tStatus\tReference\tA\n".
"chr\t2\tvalid\tA\tC\n".
"chr\t3\tvalid\tA\tC\n";

my $case2 =
"#Chromosome\tPosition\tStatus\tReference\tA\tB\n".
"chr\t2\tvalid\tA\tC\tC\n".
"chr\t3\tvalid\tA\tC\tG\n";

my $case3 =
"#Chromosome\tPosition\tStatus\tReference\tA\tB\tC\n".
"chr\t1\tvalid\tA\tA\tA\tT\n".
"chr\t2\tvalid\tA\tC\tC\tT\n".
"chr\t3\tvalid\tA\tC\tG\tT\n";

my $case4 = 
"#Chromosome\tPosition\tStatus\tA\tB\n".
"chr\t3\tvalid\tC\tG\n";

my $case5 = 
"#Chromosome\tPosition\tStatus\tB\tC\n".
"chr\t1\tvalid\tA\tT\n".
"chr\t2\tvalid\tC\tT\n".
"chr\t3\tvalid\tG\tT\n";

my $case6 = 
"#Chromosome\tPosition\tStatus\tReference\tB\n".
"chr\t2\tvalid\tA\tC\n".
"chr\t3\tvalid\tA\tG\n";

### MAIN ###

print "Testing $extract_snps_bin\n";

# writing temp file
my ($tfh,$example_file) = tempfile('extract_snps_metaalign.test.XXXXXX', TMPDIR => 1, UNLINK => $delete_temp);
print $tfh $example;
close($tfh);

run_case("Case 1", $example_file, $case1, '-r -s A');
run_case("Case 2", $example_file, $case2, '-r -s A -s B');
run_case("Case 3", $example_file, $case3, '-r -s A -s B -s C');
run_case("Case 4", $example_file, $case4, '-s A -s B');
run_case("Case 5", $example_file, $case5, '-s B -s C');
run_case("Case 6", $example_file, $case6, '-r -s B');

done_testing();
