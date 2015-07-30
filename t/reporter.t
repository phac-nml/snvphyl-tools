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

my $reporter = Reporter->new;

my $run_params = $reporter->record_run_parameters('{}', '-V', '-pe smp 4', '-pe smp 16', '-pe smp 4', '--pvar 0 --ploidy 1', '200', '15',
                                                  'mapping', '1', '-n 16 -f samsoft -r -1 -y 0.5', '-n 16 -f samsoft -r -1 -y 0.5', '-k 13 -s 6', '--pvar 0 --ploidy 1', '8',
                                                   '1', '/path/to/masked/positions.txt', 'bam1.bam', 'bam2.bam', 'bam3.bam');
is_valid_json($run_params, "The run parameters are being output in correct JSON format.");
ok((from_json($run_params))->{'parameters'}{'max_coverage'} eq '200', "The run params are correct.");
my $results = $reporter->bam_quality_data($run_params, 'reporter/sample1.bam', 'reporter/sample2.bam', 'reporter/sample3.bam');
is_valid_json($results, "The bam quality stats are being generated properly.");
ok((from_json($results))->{'bam_stats'}{'min_depth'} eq '15', "Bam stats are reported properly.");
my $filter_stats = $reporter->record_filter_stats($results, 'reporter/pseudoalign-positions1.tsv');
is_valid_json($filter_stats, "The filter stats data is being generated properly.");
ok((from_json($filter_stats))->{'filter_stats'}{'sites_unfiltered'} eq '4648', "Filter stats data is correct.");
my $file_sizes = $reporter->record_file_sizes($filter_stats, 'bam', 'reporter/sample1.bam', 'reporter/sample2.bam', 'reporter/sample3.bam');
my $file_sizes = $reporter->record_file_sizes($file_sizes, 'reference', 'reporter/reference.fasta');
is_valid_json($file_sizes, "File size data is being formatted properly.");
ok((from_json($file_sizes))->{'file_sizes'}{'bam'}{'reporter/sample1.bam'} eq '3M', "File size data is correct.");
my $ref_stats = $reporter->record_reference_info($file_sizes, 'reporter/reference.fasta', 'Illumina', 'NCBI', 'YES', 
                                                 'Escherichia', 'coli', 'O157');
is_valid_json($ref_stats, "The reference stats are being generated and formatted properly.");
ok((from_json($ref_stats))->{'reference'}{'source'} eq 'NCBI', "The reference information is correct.");
my $vcfstats = $reporter->vcf2core_stats($ref_stats, 'reporter/vcf2core.out');
is_valid_json($vcfstats, "The vcfstats are being formatted properly.");
ok((from_json($vcfstats))->{'vcf2core_stats'}{'all'}{'total_core'} eq '4910571', "The vcfstats are correct.");

#4=> test that the reference information is being recorded correctly
done_testing();
