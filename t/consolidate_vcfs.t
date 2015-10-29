#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use lib $FindBin::Bin.'/../lib';
use Test::More tests => 272;
use File::Temp 'tempfile';
use File::Temp qw /tempdir/;
use Getopt::Long;
use List::Compare;

use CompareFiles;

my $script_dir = $FindBin::Bin;

my $old_env = $ENV{'PERL5LIB'};
$ENV{'PERL5LIB'} = "$script_dir/../../../lib:$script_dir/../../../cpanlib/lib/perl5:";
$ENV{'PERL5LIB'} .= $old_env if (defined $old_env);

my $vcf_align_bin = "$script_dir/../consolidate_vcfs.pl";

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
	my ($expected,$actual) = @_;

	my $success = 1;

        #check to see if both files are empty, if so, they are the same.
        if ( (not -e $actual) && -s $expected ==0 ) {
            return $success;
        }

	### COMPARING BODY STRINGS ###
	my $expected_out = `bcftools view -H $expected`;
	my $actual_out = `bcftools view -H $actual`;

	my @expected_lines = split /\n/, $expected_out;
	my @actual_lines = split /\n/, $actual_out;

	is (scalar @actual_lines, scalar @expected_lines, 'Testing for equal body lengths.');

	my $lc = List::Compare->new(\@expected_lines, \@actual_lines);
	my @comparison_result = $lc->get_symmetric_difference();

	my $pass = 1;

	for (my $i = 0; $i <= $#comparison_result; ++$i)
	{
		my $curr_line = $comparison_result[$i];

		$pass = 0;

		if ($verbose && ($i % 2 == 0))
		{
			my $next_line = $comparison_result[$i + 1];
			print STDERR "\nMismatch in body: \n'$curr_line'\n!=\n'$next_line'\n\n";
		}
	}

	ok ($pass == 1, 'Testing bcf body strings without header.');


	### COMPARING HEADER STRINGS ###
	$expected_out = `bcftools view -h $expected`;
	$actual_out = `bcftools view -h $actual`;

	@expected_lines = split /\n/, $expected_out;
	@actual_lines = split /\n/, $actual_out;

	is (scalar @expected_lines, scalar @actual_lines, 'Testing for equal header lengths.');

	$lc = List::Compare->new(\@expected_lines, \@actual_lines);
	@comparison_result = $lc->get_symmetric_difference();

	$pass = 1;

	for (my $i = 0; $i <= $#comparison_result; ++$i)
	{
		my $curr_line = $comparison_result[$i];

		if ($curr_line !~ /##bcftools_(view|isec|filter|annotate|merge)Command=(view|isec|filter|annotate|merge .+)/)
		{
			$pass = 0;

			if ($verbose && ($i % 2 == 0))
			{
				my $next_line = $comparison_result[$i + 1];
				print STDERR "\nMismatch in headers: \n'$curr_line'\n!=\n'$next_line'\n\n";
			}
		}
	}

	ok ($pass == 1, 'Testing bcf header strings without body.');

}

sub run_command
{
	my ($freebayes,$mpileup,$coverage_cutoff,$output) = @_;

	my $command = "$vcf_align_bin --vcfsplit $freebayes --mpileup $mpileup --coverage-cutoff $coverage_cutoff --min-mean-mapping 30 --ao 0.75 --output $output";

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
}

sub get_bcfs
{
    my ($in_dir) = @_;
    my %bcf;
    my $dh;

    opendir($dh, $in_dir) or die "error opening directory $in_dir: $!";
    %bcf = map { /^(.*)\.bcf\.gz$/; $1 => "$in_dir/$_"} grep { /\.bcf\.gz$/ } readdir($dh);
    closedir($dh);

    return \%bcf
}


sub test_header
{
	my ($message) = @_;

	print "\n########################################\n";
	print "### Testing $message ###\n";
	print "########################################\n";
}




my $cases_dir = "$script_dir";
my $input_dir = "$cases_dir/consolidate_vcfs_input";

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
my @in_files = sort grep {$_ !~ /^\./} readdir($in_h);
closedir($in_h);

print "Testing all input variants in $input_dir\n\n";
for my $dir (@in_files)
{
	my $curr_input = "$input_dir/$dir";
	next if (not -d $curr_input);

	my $description = `cat $curr_input/description`;

	my $freebayes_dir = "$curr_input/input/freebayes";
	my $mpileup_dir = "$curr_input/input/mpileup";

	my $expected_output_dir = "$curr_input/output";

	my %freebayes_files = %{get_bcfs($freebayes_dir)};
	my %mpileup_files = %{get_bcfs($mpileup_dir)};

	test_header($curr_input);
	foreach my $curr_freebayes (keys %freebayes_files)
	{
		print "\n\n### Testing with $curr_freebayes.bcf.gz ###\n";
		if (exists $mpileup_files{$curr_freebayes})
		{
			my $freebayes = "$freebayes_dir/$curr_freebayes.bcf.gz";
			my $mpileup = "$mpileup_dir/$curr_freebayes.bcf.gz";
			my $expected = "$expected_output_dir/$curr_freebayes.bcf.gz";


			my $output= tempdir (CLEANUP => 1);
			$output .= "/$curr_freebayes.bcf.gz";


			run_command($freebayes,$mpileup,$coverage_cutoff,$output);



			compare_files($expected,$output);

		}
		else
		{
			die "Freebayes file $curr_freebayes.bcf.gz does not have mpileup counterpart!";
		}

	}

	my $done_testing = 0;

	print "### done ###\n";
}

done_testing();
