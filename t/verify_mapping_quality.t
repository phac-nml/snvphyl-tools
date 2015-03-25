#!/usr/bin/env perl

#Test module for the verify_mapping_quality.pl script.
use warnings;
use strict;

use Test::More;
use Test::Exception;
use File::Compare;

my ($command);
#==============================================================================
#UNIT TESTS
#=======VALID INPUT TESTS=============
#1 => verify that the script saves output in the proper file locations
system("perl ../verify_mapping_quality.pl -l mapping/output/1/ --bam bam1=mapping/input/sample1.bam --bam bam2=mapping/input/sample2.bam --bam bam3=mapping/input/sample3.bam --bam bam4=mapping/input/sample4.bam");
ok(-e 'mapping/output/1/', "Output file in correct location.");

#2 => verify that the appropriate number of isolates are logged when default values are used
system("perl ../verify_mapping_quality.pl -l mapping/output/2/ --bam bam1=mapping/input/sample1.bam --bam bam2=mapping/input/sample2.bam --bam bam3=mapping/input/sample3.bam --bam bam4=mapping/input/sample4.bam");
ok(countLines('mapping/output/2/mapping_percentage.log') == 2, "The correct number of isolates are logged when default values are used.");

#3 => verify that the appropriate number of isolates are logged when the minimum depth is changed
system("perl ../verify_mapping_quality.pl -l mapping/output/3/ --bam bam1=mapping/input/sample1.bam --bam bam2=mapping/input/sample2.bam --bam bam3=mapping/input/sample3.bam --bam bam4=mapping/input/sample4.bam --min-depth 35");
ok(countLines('mapping/output/3/mapping_percentage.log') == 3, "The correct number of isolates are logged when min_depth altered.");

#4 => verify that the appropriate number of isolates are logged when the minimum mapping percentage is changed
system("perl ../verify_mapping_quality.pl -l mapping/output/4/ --bam bam1=mapping/input/sample1.bam --bam bam2=mapping/input/sample2.bam --bam bam3=mapping/input/sample3.bam --bam bam4=mapping/input/sample4.bam -min-map 99.9");
ok(countLines('mapping/output/4/mapping_percentage.log') == 6, "The correct number of isolates are logged when min_mapping altered.");

#========INVALID INPUT TESTS========================
#5 => verify that script dies with error message when no valid bam files are input
$command = "perl ../verify_mapping_quality.pl -l mapping/output/5/ --bam nobam=mapping/no_bams_here";
ok(system(`$command`)!=0, "Invalid bam file test.");

#6 => verify that the script allows the default number of cores to be changed
$command = "perl ../verify_mapping_quality.pl -l mapping/output/6/ -i mapping/input/ --bam bam1=mapping/input/sample1.bam --bam bam2=mapping/input/sample2.bam --bam bam3=mapping/input/sample3.bam --bam bam4=mapping/input/sample4.bam -c 2";
ok(system(`$command`)!= 0, "Missing genome size test.");

#7 => verify that the script prints the log to the pwd when no logfile location is specified
system("perl ../verify_mapping_quality.pl -i mapping/input/ --bam bam1=mapping/input/sample1.bam --bam bam2=mapping/input/sample2.bam --bam bam3=mapping/input/sample3.bam --bam bam4=mapping/input/sample4.bam");
ok(-e 'mapping_percentage.log', "Log file produced in default location properly");

#Clean up any temp files produced:
system("rm mapping_percentage.log");
 
done_testing();

#method to quickly count the number of lines in output text file
sub countLines{
my ($filename) = @_;
#Count number of lines in a file
 my $lines = 0;
    open(FILE, $filename) or die "Can't open `$filename': $!";
    while (<FILE>) {
        $lines++;
    }
    close FILE;
    return $lines;
}