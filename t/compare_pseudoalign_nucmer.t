#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use File::Temp 'tempfile';
use Getopt::Long;
use Template;
use File::Basename;

my $script_dir = $FindBin::Bin;
my $compare_snps_bin = "$script_dir/../compare_snv_align_nucmer.pl";
my $delete_temp = 1;
my $verbose = 0;
my $tt = Template->new({ABSOLUTE => 1});

sub compare_files
{
        my ($expected_out_file,$actual_out_file) = @_;

        my $success = 1;

        open(my $out_h, $expected_out_file) or die "Could not open $expected_out_file: $!";
        open(my $a_out_h, $actual_out_file) or die "Could not open $actual_out_file: $!";
        while($success and (defined (my $expected_line = readline($out_h))))
        {
                my $actual_line = readline($a_out_h);
                if (not defined $actual_line)
                {
                        $success = 0;
                        fail("expected file $expected_out_file has more lines than actual file $actual_out_file");
                        next;
                }
                else
                {
                        chomp $expected_line;
                        chomp $actual_line;
                        if ($actual_line ne $expected_line)
                        {
                                is($actual_line,$expected_line,"lines \"$actual_line\" and \"$expected_line\" differ");
                                $success = 0;
                        }
                }
        }
        close($out_h);
	if (defined readline($a_out_h))
	{
		$success = 0;
		fail("actual file $actual_out_file has more lines than expected file $expected_out_file");
	}
        close($a_out_h);

        return $success;
}

sub run_case
{
	my ($name, $reference, $fasta, $snv_align, $core_positions_file, $bad_positions_file, $expected_out, $expected_detailed_out) = @_;

	my ($afh, $actual_out_file) = tempfile("compare_nucmer.actual_out.XXXXXX", TMPDIR => 1, UNLINK => $delete_temp);
	close($afh);
	print STDERR "Actual Out: $actual_out_file\n" if ($verbose);

	my ($dfh, $detailed_out_file) = tempfile("compare_nucmer.detailed_out.XXXXXX", TMPDIR => 1, UNLINK => $delete_temp);
	close($dfh);
	print STDERR "Detailed Out: $detailed_out_file\n" if ($verbose);
	
	my $command = "$compare_snps_bin -i $snv_align -r $reference -g $fasta -c $core_positions_file -b $bad_positions_file -o $detailed_out_file 1> $actual_out_file";
	$command .= " 2> /dev/null" if (not $verbose);
	print STDERR "$command\n" if ($verbose);
	system($command) == 0 or die "Could not execute command $command: $!";
	pass("pass $name") if (compare_files($expected_out, $actual_out_file) and
		compare_files($expected_detailed_out,$detailed_out_file));
}

sub build_expected_out
{
	my ($summary_template, $detailed_template, $reference, $query, $core_positions, $bad_positions, $snv_align) = @_;

	my $reference_base = basename($reference);
	my $query_base = basename($query);
	my $bad_positions_base = basename($bad_positions);
	my $core_positions_base = basename($core_positions);

	my $vars = {
		'reference' => $reference,
		'reference_base' => $reference_base,
		'query' => $query,
		'query_base' => $query_base,
		'bad_positions' => $bad_positions,
		'bad_positions_base' => $bad_positions_base,
		'core_positions' => $core_positions,
		'core_positions_base' => $core_positions_base,
		'snv_align' => $snv_align
		};

	my ($efh, $expected_summary) = tempfile("compare_snv_align.summary.XXXXXX", TMPDIR => 1, UNLINK => $delete_temp);
	close($efh);
	$tt->process($summary_template, $vars, $expected_summary) || die $tt->error(),"\n";

	my ($dfh, $expected_detailed) = tempfile("compare_nucmer.expected_detailed.XXXXXX", TMPDIR => 1, UNLINK => $delete_temp);
	close($dfh);
	$tt->process($detailed_template, $vars, $expected_detailed) || die $tt->error(),"\n";

	return ($expected_summary, $expected_detailed);
}

### MAIN ###
my ($help);
my $usage = "Usage: $0 [-h|--help] [-v|--verbose]\n";
my $tests_dir = "$script_dir/compare_snv_align_data";

if (!GetOptions('v|verbose' => \$verbose,
                'h|help' => \$help))
{
        die "Invalid option\n$usage";
}

if ($help)
{
	print $usage;
	exit 0;
}

if (not defined $verbose)
{
	$verbose = 0;
}
elsif ($verbose)
{
	$delete_temp = 0;
}

opendir(my $test_in,$tests_dir) or die "Could not open $tests_dir: $!";
my @tests = map {"$tests_dir/$_"} sort {$a <=> $b} grep {$_ !~ /^\./} readdir($test_in);
closedir($test_in);

for my $dir (@tests)
{
	print "\nTesting $dir:\n" if ($verbose);
	my $info = `cat $dir/info.txt`;
	chomp $info;

	my $reference_file = "$dir/input/reference.fasta";
	my $query_file = "$dir/input/query.fasta";
	my $snv_align_file = "$dir/input/snv_align.tsv";
	my $bad_positions_file = "$dir/input/bad_positions.tsv";
	my $core_positions_file = "$dir/input/core_positions.tsv";

	my $expected_summary_template = "$dir/output/summary.tsv.template";
	my $expected_detailed_template = "$dir/output/detailed.tsv.template";

	my ($expected_summary_file,$expected_detailed_file) = 
		build_expected_out($expected_summary_template, $expected_detailed_template,
			$reference_file, $query_file, $core_positions_file, $bad_positions_file, $snv_align_file);

	run_case("$dir: $info", $reference_file, $query_file, $snv_align_file, $core_positions_file, $bad_positions_file, $expected_summary_file,$expected_detailed_file);
}

done_testing();
