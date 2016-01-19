#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use Test::Exception;
use File::Compare;
use File::Temp qw(tempdir);

my $script_dir = $FindBin::Bin;

my $old_env = $ENV{'PERL5LIB'};
$ENV{'PERL5LIB'} = "$script_dir/../lib:$script_dir/../cpanlib/lib/perl5:";
$ENV{'PERL5LIB'} .= $old_env if (defined $old_env);

my $mapping_dir = "$script_dir/excess_coverage";
my $mapping_bin = "perl $script_dir/../verify_excess_coverage.pl";
my ($command);

#==============================================================================
#UNIT TESTS

#1 => verify that the script saves output in the proper file locations
my $create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("$mapping_bin --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam > ".$create_dir."/1.txt 2>&1");
ok(-e $create_dir."/1.txt", "Output file in correct location.");

#2 => verify that the appropriate number of isolates are logged when default values are used
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1);
system("$mapping_bin -c 8 --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam > ".$create_dir."/2.txt 2>&1");
ok(countLines($create_dir."/2.txt") == 4, "The correct number of regions are output when default values are used.");

#3 => verify that the appropriate number of isolates are logged when the max deviaiton is changed
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("$mapping_bin -c 8 --max-dev 7 --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam > ".$create_dir."/3.txt 2>&1");
ok(countLines($create_dir."/3.txt") == 5, "The correct number of isolates are logged when min_depth altered.");

#4 => verify that script dies with error message when no valid bam files are input
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
$command = "$mapping_bin -c 8 --bam bam1=$mapping_dir/input/FAKE_BAM.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam > ".$create_dir."/4.txt 2>&1";
my $return_code = system($command);
ok($return_code !=0, "Invalid bam file test.");

#5 => verify that the script allows the default number of cores to be changed
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
$command = "$mapping_bin -c 8 --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam > ".$create_dir."/5.txt 2>&1";
$return_code = system($command);
ok($return_code == 0, "Change number of cores test.");

#6 => verify that the script produces the exact output expected
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("$mapping_bin -c 8 --max-dev 3 --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam > ".$create_dir."/6.txt 2>&1");
my $strains_present = (`grep -c "#bam1" $create_dir/6.txt` && `grep -c "#bam2" $create_dir/6.txt` && `grep -c "#bam3" $create_dir/6.txt` && `grep -c "#bam4" $create_dir/6.txt`);
ok($strains_present, "The correct keys are being used to extract isolate id's.");
 
done_testing();

#method to quickly count the number of lines in output text file
sub countLines{
my ($filename) = @_;
#Count number of lines in a file
my $lines = 0;
    open(FILE, '<', $filename) or die "Can't open `$filename': $!";
    while (<FILE>) {
        $lines++;
    }
    close FILE;
    return $lines;
}