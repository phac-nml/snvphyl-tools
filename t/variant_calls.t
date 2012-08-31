#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use File::Temp 'tempfile';
use Getopt::Long;

my $script_dir = $FindBin::Bin;
my $vcf_align_bin = "$script_dir/../vcf2pseudoalignment.pl";

sub usage
{
	"Usage: $0 [Options]\n".
	"Options:\n".
	"\t-h|--help\n".
	"\t-v|--verbose\n";
}

my $cases_dir = "$script_dir";
my $input_dir = "$cases_dir/input";

my $coverage_cutoff = 4;

### MAIN ###

my ($help,$verbose);
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
my @in_files = sort grep {$_ !~ /^\./} readdir($in_h);
closedir($in_h);

print "Testing all input variants in $input_dir\n";
for my $dir (@in_files)
{
	my $curr_input = "$input_dir/$dir";
	next if (not -d $curr_input);

	my $description = `cat $curr_input/description`;
	my $expected = `cat $curr_input/expected.fasta`;
	my $reference = 'ref';
	my $expected_out_file = "$curr_input/expected.fasta";
	my $vcf_dir = $curr_input;
	my $pileup_dir = "$curr_input/pileup";
	print "\n########################################\n";
	print "### Testing $curr_input ###\n";
	print "########################################\n\n";
	print "### Description ###\n";
	print "$description\n";
	print "### Expected ###\n";
	print "$expected\n";

	die "$expected_out_file does not exist" if (not -e $expected_out_file);

	my $done_testing = 0;
	my ($fh,$actual_out_file) = tempfile('vcf2pseudoalignment.test.XXXXXXXX', TMPDIR => 1, UNLINK => 1);
	close($fh);
	my $command = "$vcf_align_bin --vcf-dir $vcf_dir --mpileup-dir $pileup_dir --reference $reference --format fasta --output $actual_out_file --coverage-cutoff $coverage_cutoff";
	
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

	my $got = `cat $actual_out_file`;
	print "### Got ###\n";
	print "$got\n";
	open(my $out_h, $expected_out_file) or die "Could not open $expected_out_file: $!";
	open(my $a_out_h, $actual_out_file) or die "Could not open $actual_out_file: $!";
	while(not $done_testing and (defined (my $expected_line = readline($out_h))))
	{
		my $actual_line = readline($a_out_h);
		if (not defined $actual_line)
		{
			$done_testing = 1;
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
				$done_testing = 1;
			}
		}
	}
	if (not $done_testing) #pass
	{
		pass("pseudoalignment generated from data in $curr_input is valid");
	}
	print "### done ###\n";
	close($out_h);
}

done_testing();
