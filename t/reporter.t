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

my $script_dir = $FindBin::Bin;

my $old_env = $ENV{'PERL5LIB'};
$ENV{'PERL5LIB'} = "$script_dir/../lib:$script_dir/../cpanlib/lib/perl5:";
$ENV{'PERL5LIB'} .= $old_env if (defined $old_env);

my $mapping_dir = "$script_dir/reporter";
my $mapping_bin = "$script_dir/../verify_mapping_quality.pl";
my $command;

#==============================================================================
#UNIT TESTS
my $create_dir = tempdir(TEMPLATE => 'tempjsonXXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";

my $output = `perl $script_dir/../reporter.pl --output=$create_dir/reporter.json --step=bam_quality_data --bam bam1=$script_dir/reporter/sample1.bam --bam bam2=$script_dir/reporter/sample2.bam --bam bam3=$script_dir/reporter/sample3.bam`;

$output = `perl $script_dir/../reporter.pl --step=record_filter_stats --output=$create_dir/reporter.json --pseudo=/t/reporter/pseudoalign-positions1.tsv --json=$create_dir/reporter.json`;
ok(check_json("$create_dir/reporter.json"), "The json for filter stats is correct.");

$output = `perl $script_dir/../reporter.pl --step=record_reference_info --output=$create_dir/reporter.json --json=$create_dir/reporter.json --ref-file='t/reporter/reference.fasta' --ref-sequencer='Illumina' --ref-source='NCBI' --plasmids='YES' --genus='Escherichia' --species='coli' --serotype='O157'`;
ok(check_json("$create_dir/reporter.json"), "The json for reference info is correct.");

$output = `perl $script_dir/../reporter.pl --step=record_file_sizes --output=$create_dir/reporter.json --json=$create_dir/reporter.json --file-type='bam' --file-sizes file1=t/reporter/sample1.bam --file-sizes file2=t/reporter/sample2.bam --file-sizes file3=t/reporter/sample3.bam`;
ok(check_json("$create_dir/reporter.json"), "The json for file sizes is correct.");

$output = `perl $script_dir/../reporter.pl --step=record_run_parameters --output=$create_dir/reporter.json --json=$create_dir/reporter.json --freebayes-params='--pvar 0 --ploidy 1' --max-coverage=200 --min-coverage=15 --processors='1' --smalt-index='-n 16 -f samsoft -r -1 -y 0.5' --smalt-map='-n 16 -f samsoft -r -1 -y 0.5' --trim-clean='-k 13 -s 6' --vcf2core-cpus=8 --run-id='123424' --masked-positions='/path/to/masked/positions.txt' --read-file read1='read1.fasta' --read-file read2='read2.fasta' --read-file read3='read3.fasta'`;
ok(check_json("$create_dir/reporter.json"), "The json for run params is correct.");

$output = `perl $script_dir/../reporter.pl --step=vcf2core_stats --output=$create_dir/reporter.json --json=$create_dir/reporter.json --vcf2core-stats=t/reporter/vcf2core.out`;
ok(check_json("$create_dir/reporter.json"), "The json for vcf2core stats is correct.");

done_testing(); 
                                                                                                    
sub check_json{
	my($json_file) = @_;
	my $json_daisy_chain="";
	open(JSON, $json_file);
	
	while(<JSON>){
   	   $json_daisy_chain .= $_;
   	}
   	
   	if(is_valid_json($json_daisy_chain, "The json string is valid.")){
   		return 1;
   	}
   	else{
   		return 0;
   	}
   	
}
