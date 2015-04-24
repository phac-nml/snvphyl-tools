#!/usr/bin/env perl

#Test module for the verify_mapping_quality.pl script.
use warnings;
use strict;

use FindBin;
use Test::More;
use Test::Exception;
use File::Compare;

my $script_dir = $FindBin::Bin;

my $old_env = $ENV{'PERL5LIB'};
$ENV{'PERL5LIB'} = "$script_dir/../lib:$script_dir/../cpanlib/lib/perl5:";
$ENV{'PERL5LIB'} .= $old_env if (defined $old_env);

my $mapping_dir = "$script_dir/mapping";
my $mapping_bin = "$script_dir/../verify_mapping_quality.pl";
my ($command);
#==============================================================================
#UNIT TESTS

#1 => verify that the script saves output in the proper file locations
system("$mapping_bin -l $mapping_dir/output/1/ --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam 2>/dev/null");
ok(-e "$mapping_dir/output/1/mapping_percentage_10x_80%.log", "Output file in correct location.");

#2 => verify that the appropriate number of isolates are logged when default values are used
system("$mapping_bin -l $mapping_dir/output/2/ --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam 2>/dev/null");
ok(countLines("$mapping_dir/output/2/mapping_percentage.log") == 2, "The correct number of isolates are logged when default values are used.");

#3 => verify that the appropriate number of isolates are logged when the minimum depth is changed
system("$mapping_bin -l $mapping_dir/output/3/ --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam --min-depth 35 2>/dev/null");
ok(countLines("$mapping_dir/output/3/mapping_percentage.log") == 6, "The correct number of isolates are logged when min_depth altered.");

#4 => verify that the appropriate number of isolates are logged when the minimum mapping percentage is changed
system("$mapping_bin -l $mapping_dir/output/4/ --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam -min-map 99.9 2>/dev/null");
ok(countLines("$mapping_dir/output/4/mapping_percentage.log") == 6, "The correct number of isolates are logged when min_mapping altered.");

#5 => verify that script dies with error message when no valid bam files are input
$command = "$mapping_bin -l $mapping_dir/output/5/ --bam nobam=$mapping_dir/no_bams_here 2>/dev/null";
my $return_code = system($command);
ok($return_code !=0, "Invalid bam file test.");

#6 => verify that the script allows the default number of cores to be changed
$command = "$mapping_bin -l $mapping_dir/output/6/ --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam -c 2 2>/dev/null";
$return_code = system($command);
ok($return_code == 0, "Change number of cores test.");

#7 => verify that the script prints the log to the pwd when no logfile location is specified
system("$mapping_bin --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam 2>/dev/null");
ok(-e "mapping_percentage_10x_80%.log", "Log file produced in default location properly");

#8 => verify that the script will work with a draft genome that contains multiple reference contigs to map against
system("$mapping_bin -l $mapping_dir/output/7/ --bam bam1=$mapping_dir/input/draft_reference/sample1.bam 2>/dev/null");
ok(-e "$mapping_dir/output/7/mapping_percentage_10x_80%.log", "Script runs when draft genomes are used.");

#9 => verify that the script return the appropriate number of samples for the draft genome
system("$mapping_bin -l $mapping_dir/output/8/ --bam bam1=$mapping_dir/input/draft_reference/sample1.bam --min-map 99.9 2>/dev/null");
ok(countLines("$mapping_dir/output/8/mapping_percentage_10x_99.9%.log") == 3, "The correct output is generated when the reads map to multiple draft reference contigs.");

#Clean up any temp files produced:
system("rm mapping_percentage_10x_80%.log");
 
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
