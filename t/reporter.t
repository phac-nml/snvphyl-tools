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
my $command;

#==============================================================================
#UNIT TESTS

my $output = `perl ../lib/reporter.pl --output=reporter.json --step=bam_quality_data --bam bam1=reporter/sample1.bam --bam bam2=reporter/sample2.bam --bam bam3=reporter/sample3.bam`;
#is_valid_json($output, "The bam quality stats are being generated properly.");

$output = `perl ../lib/reporter.pl --step=record_filter_stats --output=reporter.json --pseudo=/../t/reporter/pseudoalign-positions1.tsv --json=reporter.json`;

#is_valid_json($output, "The filter stats data is being generated properly.");

$output = `perl ../lib/reporter.pl --step=record_reference_info --output=reporter.json --json=reporter.json --ref-file='../t/reporter/reference.fasta' --ref-sequencer='Illumina' --ref-source='NCBI' --plasmids='YES' --genus='Escherichia' --species='coli' --serotype='O157'`;
#is_valid_json($output, "The reference stats are being generated and formatted properly.");

$output = `perl ../lib/reporter.pl --step=record_file_sizes --output=reporter.json --json=reporter.json --file-type='bam' --file-sizes file1=../t/reporter/sample1.bam --file-sizes file2=../t/reporter/sample2.bam --file-sizes file3=../t/reporter/sample3.bam`;
#is_valid_json($output, "File size data is being formatted properly.");

$output = `perl ../lib/reporter.pl --step=record_run_parameters --output=reporter.json --json=reporter.json --freebayes-params='--pvar 0 --ploidy 1' --max-coverage=200 --min-coverage=15 --mode='mapping' --processors='1' --smalt-index='-n 16 -f samsoft -r -1 -y 0.5' --smalt-map='-n 16 -f samsoft -r -1 -y 0.5' --trim-clean='-k 13 -s 6' --vcf2core-cpus=8 --run-id='123424' --masked-positions='/path/to/masked/positions.txt' --read-file read1='read1.fasta' --read-file read2='read2.fasta' --read-file read3='read3.fasta'`;
#is_valid_json($output, "The run parameters are being output in correct JSON format.");

$output = `perl ../lib/reporter.pl --step=vcf2core_stats --output=reporter.json --json=reporter.json --vcf2core-stats=../t/reporter/vcf2core.out`;
#is_valid_json($output, "vcf2core_stats json is being generated correctly.");

done_testing(); 

                                                                                                      
#is_valid_json($run_params, "The run parameters are being output in correct JSON format.");
#ok((from_json($run_params))->{'parameters'}{'max_coverage'} eq '200', "The run params are correct.");
#my $results = $reporter->bam_quality_data($run_params, 'reporter/sample1.bam', 'reporter/sample2.bam', 'reporter/sample3.bam');
#is_valid_json($results, "The bam quality stats are being generated properly.");
#ok((from_json($results))->{'bam_stats'}{'min_depth'} eq '15', "Bam stats are reported properly.");
#my $filter_stats = $reporter->record_filter_stats($results, 'reporter/pseudoalign-positions1.tsv');
#is_valid_json($filter_stats, "The filter stats data is being generated properly.");
#ok((from_json($filter_stats))->{'filter_stats'}{'sites_unfiltered'} eq '4648', "Filter stats data is correct.");
#my $file_sizes = $reporter->record_file_sizes($filter_stats, 'bam', 'reporter/sample1.bam', 'reporter/sample2.bam', 'reporter/sample3.bam');
#my $file_sizes = $reporter->record_file_sizes($file_sizes, 'reference', 'reporter/reference.fasta');
#is_valid_json($file_sizes, "File size data is being formatted properly.");
#ok((from_json($file_sizes))->{'file_sizes'}{'bam'}{'reporter/sample1.bam'} eq '3M', "File size data is correct.");
#my $ref_stats = $reporter->record_reference_info($file_sizes, 'reporter/reference.fasta', 'Illumina', 'NCBI', 'YES', 
#                                                 'Escherichia', 'coli', 'O157');
#is_valid_json($ref_stats, "The reference stats are being generated and formatted properly.");
#ok((from_json($ref_stats))->{'reference'}{'source'} eq 'NCBI', "The reference information is correct.");
#my $vcfstats = $reporter->vcf2core_stats($ref_stats, 'reporter/vcf2core.out');
#is_valid_json($vcfstats, "The vcfstats are being formatted properly.");
#ok((from_json($vcfstats))->{'vcf2core_stats'}{'all'}{'total_core'} eq '4910571', "The vcfstats are correct.");

#4=> test that the reference information is being recorded correctly

