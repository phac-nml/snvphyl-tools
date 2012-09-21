#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use File::Temp 'tempfile';
use Getopt::Long;

my $script_dir = $FindBin::Bin;
my $vcf_align_bin = "$script_dir/../vcf2pseudoalignment.pl";

my $verbose = 0;

sub usage
{
	"Usage: $0 [Options]\n".
	"Options:\n".
	"\t-h|--help\n".
	"\t-v|--verbose\n";
}

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
	close($a_out_h);

	return $success;
}

sub run_command
{
	my ($vcf_dir,$pileup_dir,$reference,$coverage_cutoff,$format) = @_;

	my ($fh,$actual_out_file) = tempfile('vcf2pseudoalignment.test.XXXXXXXX', TMPDIR => 1, UNLINK => 1);
	close($fh);
	my $command = "$vcf_align_bin --vcf-dir $vcf_dir --mpileup-dir $pileup_dir --reference $reference --format $format --output $actual_out_file --coverage-cutoff $coverage_cutoff";
	
	if ($verbose)
	{
		$command .= " -v";
	}
	else
	{
		$command .= " 2>&1 1>/dev/null";
	}
	print "## Running $command\n\n";
	system($command) == 0 or die "Could not run command $command: $!";

	return $actual_out_file;
}

sub test_header
{
	my ($message) = @_;

	print "\n########################################\n";
	print "### Testing $message ###\n";
	print "########################################\n\n";
}

my $cases_dir = "$script_dir";
my $input_dir = "$cases_dir/input";

my $coverage_cutoff = 4;

### MAIN ###

my ($help);
if (!GetOptions('help|h' => \$help,
                'verbose|v' => \$verbose))
{
        die "Invalid option\n".usage;
}

if ($help)
{
	print usage;
	exit 0;
}

opendir(my $in_h,$input_dir) or die "Could not open $input_dir: $!";
my @in_files = sort {$a <=> $b} grep {$_ !~ /^\./} readdir($in_h);
closedir($in_h);

print "Testing all input variants in $input_dir\n";
for my $dir (@in_files)
{
	my $curr_input = "$input_dir/$dir";
	next if (not -d $curr_input);

	my $description = `cat $curr_input/description`;
	my $expected = `cat $curr_input/expected.fasta`;
	my $expected_out_file = "$curr_input/expected.fasta";
	my $vcf_dir = $curr_input;
	my $pileup_dir = "$curr_input/pileup";
	test_header($curr_input);
	print "### Description ###\n";
	print "$description\n";
	print "### Expected ###\n";
	print "$expected\n";

	die "$expected_out_file does not exist" if (not -e $expected_out_file);

	my $done_testing = 0;

	my $actual_out_file = run_command($vcf_dir,$pileup_dir,'ref',$coverage_cutoff,'fasta');
	my $got = `cat $actual_out_file`;
	print "### Got ###\n";
	print "$got\n";
	my $success = compare_files($expected_out_file,$actual_out_file);
	if ($success) #pass
	{
		pass("pseudoalignment generated from data in $curr_input is valid");
	}
	print "### done ###\n";
}

my $actual_file;
my $expected_file;
my $curr_input;
my $got;
my $expected;

$curr_input = "$input_dir/1";
test_header("phylip output format in $curr_input");
$expected_file = "$curr_input/expected.phy";
$expected = `cat $expected_file`;
print "### Expected ###\n";
print "$expected\n";
die("could not find input dir $curr_input") if (not -e $curr_input);
$actual_file = run_command($curr_input,"$curr_input/pileup",'ref',$coverage_cutoff,'phylip');
$got = `cat $actual_file`;
print "### Got ###\n";
print "$got\n";
pass("pass test for phylip output") if (compare_files($expected_file,$actual_file));

done_testing();
