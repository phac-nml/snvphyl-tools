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

#1 => verify that the plugin is creating the expected output
my $create_dir = tempdir(TEMPLATE => 'tempXXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
my $regions_output = "$create_dir/density_regions.tsv";
system("bcftools plugin filter_snp_density $density_dir/input/1.vcf -O b -o $create_dir/temp1.bcf -- --filename $density_dir/input/1.vcf --region_file $regions_output  --threshold 10");
ok(-e $create_dir."/temp1.bcf", "Output bcf file is in the correct location.");
ok(-e "$create_dir/density_regions.tsv", "The output density regions txt file is being created." );
#2 => verify that negative threshold values are set to the default threshold value
system("bcftools plugin filter_snp_density $density_dir/input/1.vcf -O b -o $create_dir/temp2.bcf -- -f $density_dir/input/1.vcf --region_file $regions_output --threshold 5");
my $lines1 = `bcftools view $density_dir/input/1.vcf | wc -l 2>&1`;
my $lines2 = `bcftools view $create_dir/temp2.bcf | wc -l 2>&1`;
my $output = `bcftools view $create_dir/temp2.bcf`;
ok( $lines1 == ($lines2 - 3), "The correct number of isolates are logged when default values are used.");
ok($output =~ 'filtered-density', "SNP's are correctly being filtered for density.");
ok($output =~ '##FILTER=<ID=filtered-density,Description="Set true if spacing is < 5 bp">', "Filter descripition is added to header.");

ok(`bcftools view $create_dir/temp2.bcf | grep 'filtered-density' | wc -l` == 15, "The correct number of density positions are being filtered.");
ok(`grep '20\t80' $create_dir/density_regions.tsv | wc -l` > 0, "The bcftools plugin is assigning the correct chromosome id to positions.");
done_testing();
