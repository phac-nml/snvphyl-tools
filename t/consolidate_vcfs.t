#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use lib $FindBin::Bin.'/../lib';
use Test::More tests => 278;
use File::Temp 'tempfile';
use File::Temp qw /tempdir/;
use Getopt::Long;
use Text::Diff;
use File::Basename;

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

sub compare_positions
{
	my ($expected,$actual) = @_;
	my $output; # string holding the diff output
	my @output; # array holding the split output string (split by \n)
	my @filtered_output; # array containing @output but without the diff fluff

	open (FH, "< $expected") or die "Can't open $expected for read: $!";
	my @expected = <FH>;
	close FH or die "Cannot close $expected: $!";

	open (XH, "< $actual") or die "Can't open $actual for read: $!";
	my @actual = <XH>;
	close XH or die "Cannot close $actual: $!";

	is (scalar @actual, scalar @expected, 'Testing for equal size density-filtered positions files.');

	diff \@expected, \@actual, {OUTPUT => \$output, STYLE => "OldStyle"};

	if ($output)
	{
		@output = split /\n/, $output;
	}
	is(scalar @output, 0, 'Testing that actual density-filtered position file has no diff from expected.');

	# if differences were found, list them
	if ($verbose)
	{
		if (@output)
		{
			## Filtering output (removing diff formatting fluff)
			@filtered_output = ();
			foreach my $curr_line (@output)
			{
				if ($curr_line =~ /ref/)
				{
					push (@filtered_output, substr $curr_line, 2);
				}
			}

			foreach my $curr_line (@filtered_output)
			{
					print STDERR "Mismatch!\n'$curr_line' absent from expected output!\n\n";
			}
		}
	}
}

sub compare_bcfs
{
	my ($expected,$actual) = @_;

	my $success = 1; # final result of comparison
	my $output; # string holding the diff output
	my @output; # array holding the split output string (split by \n)
	my @filtered_output; # array containing @output but without the diff fluff

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

	diff \$expected_out, \$actual_out, {OUTPUT => \$output, STYLE => "OldStyle"};

	if ($output)
	{
		@output = split /\n/, $output;
	}
	is(scalar @output, 0, 'Testing bcf body strings without header.');

	# if differences were found, list them
	if ($verbose)
	{
		if (@output)
		{
			## Filtering output (removing diff formatting fluff)
			@filtered_output = ();
			foreach my $curr_line (@output)	{
				if ($curr_line =~ /ref/)
				{
					push (@filtered_output, substr $curr_line, 2);
				}
			}

			foreach my $curr_line (@filtered_output)
			{
				# A limitation of this implementation is that if we get mismatched command type and command, it
				# will not get caught as an error. For example: ##bcftools_viewCommand=Merge
				if ($curr_line !~ /##bcftools_(view|isec|filter|annotate|merge|plugin)Command=(view|isec|filter|annotate|merge|plugin .+)/)
				{
					# Since @comparison_result is going to be a list of pairs of differences, since we want to
					# compare the first of the pair with the second of the pair in the output, we need to make sure
					# that obj 1 is in the same pair as obj 2. We use the parity of the index to determine whether or not
					# we are in a position to compare i and i + 1.
					print STDERR "Mismatch!\n'$curr_line' absent from expected output!\n\n";
				}
			}
		}
	}



	### COMPARING HEADER STRINGS ###
	$expected_out = `bcftools view -h $expected`;
	$actual_out = `bcftools view -h $actual`;

	@expected_lines = split /\n/, $expected_out;
	@actual_lines = split /\n/, $actual_out;

	is (scalar @expected_lines, scalar @actual_lines, 'Testing for equal header lengths.');

	diff \$expected_out, \$actual_out, {OUTPUT => \$output, STYLE => "OldStyle"};
	@output = split /\n/, $output;

	## Filtering output (removing diff formatting fluff)
	@filtered_output = ();
	foreach my $curr_line (@output)	{
		if ($curr_line =~ /##/)
		{
			push (@filtered_output, substr $curr_line, 2);
		}
	}

	my $pass = 1;

	foreach my $curr_line (@filtered_output)
	{
		# A limitation of this implementation is that if we get mismatched command type and command, it
		# will not get caught as an error. For example: ##bcftools_viewCommand=Merge
		if ($curr_line !~ /##bcftools_(view|isec|filter|annotate|merge|plugin)(Version|Command)=(1.3.+|view|isec|filter|annotate|merge|plugin .+)/ && $curr_line !~ /##FILTER=/)
		{
			$pass = 0;

			# Since @comparison_result is going to be a list of pairs of differences, since we want to
			# compare the first of the pair with the second of the pair in the output, we need to make sure
			# that obj 1 is in the same pair as obj 2. We use the parity of the index to determine whether or not
			# we are in a position to compare i and i + 1.
			if ($verbose)
			{
				print STDERR "Mismatch!\n'$curr_line' absent from expected output!\n\n";
			}
		}
	}

	ok ($pass == 1, 'Testing bcf header strings without body.');
}

sub run_command
{
	my ($freebayes,$mpileup,$coverage_cutoff,$output,$test_type) = @_;

	my $filtered_density_out = dirname($output)."/density_filtered_positions.tsv";
	my $command;

	if ($test_type == 1) #simple test
	{
		$command = "$vcf_align_bin --vcfsplit $freebayes --mpileup $mpileup --coverage-cutoff $coverage_cutoff --min-mean-mapping 30 --snv-abundance-ratio 0.75 --output $output --filtered-density-out $filtered_density_out --window-size 500 --density-threshold 2";
	}
	if ($test_type == 2) #snv_density-specific test
	{
		$command = "$vcf_align_bin --vcfsplit $freebayes --mpileup $mpileup --coverage-cutoff $coverage_cutoff --min-mean-mapping 30 --snv-abundance-ratio 0.75 --output $output --filtered-density-out $filtered_density_out --window-size 200 --density-threshold 4";
	}

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
my $regular_test = 1;
my $snv_test = 2;

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
	next if ($dir eq "snv-specific" || not -d $curr_input);
	my $basename = basename($curr_input);

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


			run_command($freebayes,$mpileup,$coverage_cutoff,$output, $regular_test);

			compare_bcfs($expected,$output);
		}
		else
		{
			die "Freebayes file $curr_freebayes.bcf.gz does not have mpileup counterpart!";
		}

	}

	print "### done ###\n";
}

$input_dir = "$cases_dir/consolidate_vcfs_input/snv-specific";

opendir($in_h,$input_dir) or die "Could not open $input_dir: $!";
@in_files = sort grep {$_ !~ /^\./} readdir($in_h);
closedir($in_h);

print "Testing all input variants in $input_dir\n\n";
for my $dir (@in_files)
{
	my $curr_input = "$input_dir/$dir";
	next if (not -d $curr_input);

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
			my $expected = "$expected_output_dir/expected_density_filtered_positions.tsv";

			my $output= tempdir (CLEANUP => 1);
			$output .= "/$curr_freebayes.bcf.gz";

			run_command($freebayes,$mpileup,$coverage_cutoff,$output, $snv_test);

			my $actual = dirname($output)."/density_filtered_positions.tsv";

			compare_positions($expected,$actual);
		}
		else
		{
			die "Freebayes file $curr_freebayes.bcf.gz does not have mpileup counterpart!";
		}

	}

	print "### done ###\n";
}

done_testing();
