#!/usr/bin/env perl

use Test::More;
use Test::Exception;
use File::Compare;
use File::Temp qw(tempdir);
use FindBin;
use lib $FindBin::Bin.'/../lib';

my $create_dir = tempdir(TEMPLATE => 'tempjsonXXXX', CLEANUP => 1) or die "Unable to create a temporary file directory.";

#esnure that the results are consistent
my $results1 = `perl ../filter-stats.pl -a -i filter-stats/input1.tsv > $create_dir/temp1.txt`; 
my $diff_results = `diff filter-stats/expected1.txt $create_dir/temp1.txt`;
ok(!$diff_results);
ok(`grep "Number of sites used to generate phylogeny: 8934" $create_dir/temp1.txt`);
ok(`grep "Total number of sites identified: 9623" $create_dir/temp1.txt`);
ok(`grep "Percentage of sites filtered: 7.16" $create_dir/temp1.txt`);
ok(`grep "Coverage filtered: 0" $create_dir/temp1.txt`);
ok(`grep "mpileup filtered: 0" $create_dir/temp1.txt`);
ok(`grep "Number of sites filtered: 689" $create_dir/temp1.txt`);



my $results2 =  `perl ../filter-stats.pl -i filter-stats/input1.tsv > $create_dir/temp2.txt`; 
my $diff_without_a = `diff filter-stats/expected1.txt $create_dir/temp2.txt`;
ok($diff_without_a);

#Tests to ensure that the -a flag functions properly
my $grep_invalids = `grep "Invalid filtered: 689" $create_dir/temp1.txt`;
my $grep_invalids2 = `grep "Invalid filtered: Invalid positions not analyzed." $create_dir/temp2.txt`;
ok($grep_invalids);
ok($grep_invalids2);

 
