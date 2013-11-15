#!/usr/bin/env perl

use warnings;
use strict;

use Test::More;
use Test::Deep;
use File::Temp 'tempfile';
use YAML::Tiny;

use FindBin;
use lib $FindBin::Bin.'/../lib';

use Align::Nucmer;

my $script_dir = $FindBin::Bin;
my $verbose = 0;

### MAIN ###
my $tests_dir = "$script_dir/nucmer_align_data";

opendir(my $td,$tests_dir) or die "Could not open $tests_dir: $!";
my @tests = map {"$tests_dir/$_"} sort {$a <=> $b} grep {$_ !~ /^\./} readdir($td);
closedir($td);

my $nucmer_align_parser = Align::Nucmer->new;

for my $dir (@tests)
{
	print "\nTesting $dir:\n" if ($verbose);
	my $info_file = "$dir/info.txt";
	my $reference_file = "$dir/reference.fasta";
	my $query_file = "$dir/query.fasta";
	my $expected_file = "$dir/expected.yaml";

	print `cat $info_file` if ($verbose);

	my $yaml = YAML::Tiny->read($expected_file);
	die "error: yaml file $expected_file not valid" if (not defined $yaml);
	my $expected_result = $yaml->[0];

	my $actual_result = $nucmer_align_parser->align_and_parse($reference_file,$query_file);

	cmp_deeply($actual_result,$expected_result, "data structures okay for $dir");
}

done_testing();
