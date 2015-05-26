#!/usr/bin/env perl

#Test module for the verify_mapping_quality.pl script.
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

my $mapping_dir = "$script_dir/mapping";
my $mapping_bin = "$script_dir/../verify_mapping_quality.pl";
my ($command);

#==============================================================================
#UNIT TESTS

#1 => verify that the script saves output in the proper file locations
my $create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("$mapping_bin -l ".$create_dir."/ --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam 2>/dev/null");
ok(-e $create_dir."/mapping_percentage.log", "Output file in correct location.");

#2 => verify that the appropriate number of isolates are logged when default values are used
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1);
system("$mapping_bin -l ".$create_dir."/ --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam 2>/dev/null");
ok(countLines($create_dir."/mapping_percentage.log") == 4, "The correct number of isolates are logged when default values are used.");

#3 => verify that the appropriate number of isolates are logged when the minimum depth is changed
$create_dir = File::Temp->newdir(TEMPLATE => 'XXXXX');
system("$mapping_bin -l ".$create_dir."/ --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam --min-depth 35 2>/dev/null");
ok(countLines($create_dir->dirname."/mapping_percentage.log") == 8, "The correct number of isolates are logged when min_depth altered.");

#4 => verify that the appropriate number of isolates are logged when the minimum mapping percentage is changed
$create_dir = File::Temp->newdir(TEMPLATE => 'XXXXX');
system("$mapping_bin -l ".$create_dir."/ --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam -min-map 99.9 2>/dev/null");
ok(countLines($create_dir->dirname."/mapping_percentage.log") == 8, "The correct number of isolates are logged when min_mapping altered.");

#5 => verify that script dies with error message when no valid bam files are input
$create_dir = File::Temp->newdir(TEMPLATE => 'XXXXX');
$command = "$mapping_bin -l ".$create_dir."/ --bam nobam=$mapping_dir/no_bams_here 2>/dev/null";
my $return_code = system($command);
ok($return_code !=0, "Invalid bam file test.");

#6 => verify that the script allows the default number of cores to be changed
$create_dir = File::Temp->newdir(TEMPLATE => 'XXXXX');
$command = "$mapping_bin -l ".$create_dir."/ --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam -c 2 2>/dev/null";
$return_code = system($command);
ok($return_code == 0, "Change number of cores test.");

#7 => verify that the script prints the log to the pwd when no logfile location is specified
$create_dir = File::Temp->newdir(TEMPLATE => 'XXXXX');
system("$mapping_bin --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam 2>/dev/null");
ok(-e "mapping_percentage.log", "Log file produced in default location properly");

#8 => verify that the script will work with a draft genome that contains multiple reference contigs to map against
$create_dir = File::Temp->newdir(TEMPLATE => 'XXXXX');
system("$mapping_bin -l ".$create_dir."/ --bam bam1=$mapping_dir/input/draft_reference/sample1.bam 2>/dev/null");
ok(-e $create_dir->dirname."/mapping_percentage.log", "Script runs when draft genomes are used.");

#9 => verify that the script return the appropriate number of samples for the draft genome
$create_dir = File::Temp->newdir(TEMPLATE => 'XXXXX');
system("$mapping_bin -l ".$create_dir."/ --bam bam1=$mapping_dir/input/draft_reference/sample1.bam --min-map 99.9 2>/dev/null");
ok(countLines($create_dir->dirname."/mapping_percentage.log") == 5, "The correct output is generated when the reads map to multiple draft reference contigs.");

#Clean up any temp files produced and remove the output folder:
system("rm mapping_percentage.log");
 
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
