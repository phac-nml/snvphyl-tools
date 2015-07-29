#!/usr/bin/env perl

use FindBin;
use Test::More;
use Test::Exception;
use File::Compare;
use File::Temp qw(tempdir);
use JSON;
use Test::JSON;
use FindBin;
use lib $FindBin::Bin.'/../lib';
use Reporter;

my $script_dir = $FindBin::Bin;

my $old_env = $ENV{'PERL5LIB'};
$ENV{'PERL5LIB'} = "$script_dir/../lib:$script_dir/../cpanlib/lib/perl5:";
$ENV{'PERL5LIB'} .= $old_env if (defined $old_env);

my $mapping_dir = "$script_dir/reporter";
my $mapping_bin = "$script_dir/../verify_mapping_quality.pl";
my ($command);

#==============================================================================
#UNIT TESTS

my %test_json;

#1 => test the record_read_mapping submodule:
my $reporter = Reporter->new;
my $run_params = $reporter->record_run_parameters('{}', '-V', '-pe smp 16', '-pe smp 16', '-pe smp 4', '--pvar 0 --ploidy 1', '200', '15', 'mapping', '1', '--pvar 0 --ploidy 1', '--pvar 0 --ploidy 1', '--pvar 0 --ploidy 1', '8', '1');
my $results = $reporter->bam_quality_data($run_params, 'reporter/sample1.bam', 'reporter/sample2.bam', 'reporter/sample3.bam');

#2=> test the record_reference_info submodule:
#my $reference = $reporter->record_reference_info(to_json($test_json), '12345', 'Escherichia', 'coli', 'O157', 'Illumina', 'NML', 'NO', '3500000', 'NO');

my $filter_stats = $reporter->record_filter_stats($results, 'reporter/pseudoalign-positions1.tsv');

my $file_sizes = $reporter->record_file_sizes('bam', $filter_stats, 'reporter/sample1.bam', 'reporter/sample2.bam', 'reporter/sample3.bam');
my $file_sizes = $reporter->record_file_sizes('reference', $file_sizes, '../t/reporter/reference.fasta');
my $ref_stats = $reporter->record_reference_info($file_sizes, '../t/reporter/reference.fasta', 'Illumina', 'NCBI', 'YES');
my $vcfstats = $reporter->vcf2core_stats($ref_stats, 'reporter/vcf2core.out');

print $vcfstats;

#4=> test that the reference information is being recorded correctly

done_testing();
