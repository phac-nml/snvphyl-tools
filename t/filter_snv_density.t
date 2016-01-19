#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use Test::Exception;
use File::Compare;
use File::Temp qw(tempdir);

my $script_dir = $FindBin::Bin;
my $density_dir = "$script_dir/snv_density";
my ($command);
#=============================================================================
#UNIT TESTS

#1 => verify that the plugin is creating the expected output
my $create_dir = tempdir(TEMPLATE => 'tempXXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
my $regions_output = "$create_dir/density_regions.tsv";
my $regions_output2 = "$create_dir/density_regions2.tsv";
my $regions_output3 = "$create_dir/density_regions3.tsv";

my $result = `bcftools plugin filter_snv_density $density_dir/input/1.vcf -O b -o $create_dir/temp1.bcf -- --filename $density_dir/input/1.vcf --region_file $regions_output --window_size 150 --threshold 2`;
ok(-e $create_dir."/temp1.bcf", "Output bcf file is in the correct location.");
ok(-e $regions_output, "The output density regions txt file is being created." );
ok((`grep '2\t101\t115' $regions_output | wc -l` > 0) && (`grep '20\t3\t95' $regions_output | wc -l` > 0) &&  (`grep '20\t273\t278' $regions_output | wc -l` > 0), "The bcftools plugin is assigning the correct chromosome id to positions.");


$result = `bcftools plugin filter_snv_density $density_dir/input/1.vcf -O b -o $create_dir/temp2.bcf -- --filename $density_dir/input/1.vcf --region_file $regions_output2 --window_size 150 --threshold 4`;
ok((`grep '2\t101\t115' $regions_output2 | wc -l` == 0) && (`grep '20\t3\t95' $regions_output2 | wc -l` > 0) &&  (`grep '20\t273\t278' $regions_output2 | wc -l` == 0), "Threshold edge case works properly.");

$result = `bcftools plugin filter_snv_density $density_dir/input/1.vcf -O b -o $create_dir/temp2.bcf -- --filename $density_dir/input/1.vcf --region_file $regions_output3 --window_size 300  --threshold 2`;
ok((`grep '2\t101\t115' $regions_output3 | wc -l` == 1) && (`grep '20\t3\t278' $regions_output3 | wc -l` == 1) &&  (`grep '3\t10\t15' $regions_output3 | wc -l` == 1), "Window edge case works properly.");

done_testing();
