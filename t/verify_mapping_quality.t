#!/usr/bin/env perl

#Test module for the verify_mapping_quality.pl script.
use warnings;
use strict;

use FindBin;
use Test::More;
use Test::Exception;
use File::Compare;
use File::Temp qw(tempfile);

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
my $result = `$mapping_bin  --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam`;
ok($result, "Output is generated.");

#2 => verify that the appropriate number of isolates are logged when default values are used
$result = `$mapping_bin --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam`;
ok($result =~ tr/\n// == 4, "The correct number of isolates are logged when default values are used.");

#3 => verify that the appropriate number of isolates are logged when the minimum depth is changed
$result = `$mapping_bin --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam --min-depth 35`;
ok($result =~ tr/\n// == 7, "The correct number of isolates are logged when min_depth altered.");

#4 => verify that the appropriate number of isolates are logged when the minimum mapping percentage is changed
$result = `$mapping_bin --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam -min-map 99.9`;
ok($result =~ tr/\n// == 8, "The correct number of isolates are logged when min_mapping altered.");

#5 => verify that script dies with error message when no valid bam files are input
$command = system("$mapping_bin --bam nobam=$mapping_dir/no_bams_here 2>&1");
my $return_code = system($command);
ok($return_code !=0, "Invalid bam file test.");

#6 => verify that the script allows the default number of cores to be changed
$command = "$mapping_bin --bam bam1=$mapping_dir/input/sample1.bam --bam bam2=$mapping_dir/input/sample2.bam --bam bam3=$mapping_dir/input/sample3.bam --bam bam4=$mapping_dir/input/sample4.bam -c 2 2>&1";
$return_code = system($command);
ok($return_code == 0, "Change number of cores test.");

#7 => verify that the script will work with a draft genome that contains multiple reference contigs to map against
$result = `$mapping_bin --bam bam1=$mapping_dir/input/draft_reference/sample1.bam 2>&1`;
ok($result, "Script runs when draft genomes are used.");

#8 => Using designed output, ensure that with min-mpa at 80 and min depth at 2, 7.29% of the positions map
$result = `$mapping_bin --bam bam1=$mapping_dir/input/sample5.bam --min-depth 2 --min-map 80 2>&1`;
ok($result =~ '7.29%', "The script generates the correct mapping percentage with designed data.");


#9 => Using designed output, ensure that with min-mpa at 80 and min depth at 2, 7.29% of the positions map and using ---output option
my (undef,$output) = tempfile();
my $expected_out = "$mapping_dir/output/output_9.txt";

$command = "$mapping_bin --bam bam1=$mapping_dir/input/sample5.bam --min-depth 2 --min-map 80 --output $output 2>&1";
$return_code = system($command);
ok($return_code == 0, "Script ran successfully");
ok(compare($output,$expected_out) == 0, "Bam1 had 7.29% core with 2x on 80% of the genome");
done_testing();

