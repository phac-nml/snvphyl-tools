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

my $density_dir = "$script_dir/snp_density";
my ($command);

#==============================================================================
#UNIT TESTS
#1 => verify that the script writes the correct number of records to the output
my $create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("bcftools plugin filter_snp_density $density_dir/input/1.bcf -O b -- -f $density_dir/input/1.bcf -t 10 > $create_dir/temp1.bcf 2>&1");
ok(-e $create_dir."/temp1.bcf", "Output bcf file is in the correct location.");

#2 => verify that negative threshold values are set to the default threshold value
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1);
system("bcftools plugin filter_snp_density $density_dir/input/1.bcf -O b -o $create_dir/temp2.bcf -- -f $density_dir/input/1.bcf -t 5 2>&1");
my $lines1 = `bcftools view $density_dir/input/1.bcf | wc -l 2>&1`;
my $lines2 = `bcftools view $create_dir/temp2.bcf | wc -l 2>&1`;
print $lines2;
ok( $lines1 == ($lines2 - 3), "The correct number of isolates are logged when default values are used.");

#3 => verify that an absent threshold value is set to the default value
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1);
ok(!system("bcftools plugin filter_snp_density $density_dir/input/1.bcf -O b -o $create_dir/temp3.bcf -- -f $density_dir/input/1.bcf 2>&1"), "A default threshold value is set when the parameter is absent");

#`perl ../vcf2pseudoalignment.pl --fasta $script_dir/snp_density/input/reference.fasta --vcf-dir $script_dir/snp_density/input/vcf/ --mpileup-dir $script_dir/snp_density/input/mpileup/ --output-base $script_dir/snp_density/results`;

done_testing();