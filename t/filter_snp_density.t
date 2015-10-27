#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;
use Test::More;
use Test::Exception;
use File::Compare;
use File::Temp qw(tempdir);

my $script_dir = $FindBin::Bin;
my $density_dir = "$script_dir/snp_density";
my $regions_output = "$density_dir/density_regions.tsv";
my ($command);
#=============================================================================
#UNIT TESTS
my $create_dir = tempdir(TEMPLATE => 'tempXXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("bcftools plugin filter_snp_density $density_dir/input/1.bcf -O b -- --filename $density_dir/input/1.bcf --region_file $regions_output  --threshold 10 > $density_dir/temp1.bcf 2>&1");
ok(-e $create_dir."/temp1.bcf", "Output bcf file is in the correct location.");


#2 => verify that negative threshold values are set to the default threshold value
system("bcftools plugin filter_snp_density $density_dir/input/1.bcf -O b -o $create_dir/temp2.bcf -- -f $density_dir/input/1.bcf --region_file $regions_output --threshold 5 2>&1");
my $lines1 = `bcftools view $density_dir/input/1.bcf | wc -l 2>&1`;
my $lines2 = `bcftools view $create_dir/temp2.bcf | wc -l 2>&1`;
my $output = `bcftools view $create_dir/temp2.bcf`;
ok( $lines1 == ($lines2 - 3), "The correct number of isolates are logged when default values are used.");
ok($output =~ 'filtered-density', "SNP's are correctly being filtered for density.");
ok($output =~ '##FILTER=<ID=filtered-density,Description="Set true if spacing is < 5 bp">', "Filter descripition is added to header.");

#3 => verify that an absent threshold value is set to the default value
ok(!system("bcftools plugin filter_snp_density $density_dir/input/1.bcf -O b -o $create_dir/temp3.bcf -- -f $density_dir/input/1.bcf --region_file $regions_output 2>&1"), "A default threshold value is set when the parameter is absent");


ok(system("perl $script_dir/../filter_snp_density.pl --output-file $density_dir/new_out.bcf --output-regions $density_dir/new_regions.tsv --threshold 10 --bcf bcf1=$density_dir/input/1.bcf 2>&1"));


done_testing();
