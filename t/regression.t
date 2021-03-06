#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use lib $FindBin::Bin.'/../lib';
use Test::More;
use File::Temp 'tempfile'; 
use Getopt::Long;

use CompareFiles;

my $script_dir = $FindBin::Bin;

my $old_env = $ENV{'PERL5LIB'};
$ENV{'PERL5LIB'} = "$script_dir/../../../lib:$script_dir/../../../cpanlib/lib/perl5:";
$ENV{'PERL5LIB'} .= $old_env if (defined $old_env);

my $vcf_align_bin = "$script_dir/../bin/vcf2snv_alignment.pl";

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

        #check to see if both files are empty, if so, they are the same.
        if ( (not -e $actual_out_file) and (not -e $expected_out_file) ) {
            return $success;
        }

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
		fail("expected file $expected_out_file has less lines than actual file $actual_out_file");
	}

	close($a_out_h);

	unlink($actual_out_file);

	return $success;
}

sub run_command
{
	my ($vcf_hash, $reference, $formats, $extra_params) = @_;

	my %vcf_files = %$vcf_hash;
        
	my ($fh,$actual_out_base) = tempfile('vcf2snv_alignment.test.XXXXXXXX', TMPDIR => 1, UNLINK => 1);

        close($fh);
	my $format = '';
	$extra_params = '' if (not defined $extra_params);
	my @out_files = ();
	for my $f (@$formats)
	{
		$format .= "--format $f ";
		if ($f eq 'phylip')
		{
			push(@out_files,"$actual_out_base.phy");
		}
		elsif ($f eq 'fasta')
		{
			push(@out_files, "$actual_out_base.fasta");
		}
		else
		{
			die "Invalid format $f for testing";
		}
	}

	my $vcf_cmd = "--consolidate_vcf" . join (" --consolidate_vcf" , map { " $_=" . $vcf_files{$_} } keys %vcf_files);


	my $command = "$vcf_align_bin $vcf_cmd --reference $reference $format --output-base $actual_out_base $extra_params";


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
	return ($actual_out_base,@out_files);
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
	print "########################################\n\n";
}

my $cases_dir = "$script_dir";
my $input_dir = "$cases_dir/regression";

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

print "Testing all input variants in $input_dir\n";
for my $dir (@in_files)
{
	my $curr_input = "$input_dir/$dir";
	# curr_input = snvphyl-tools/t/input/1,2,etc.
	next if (not -d $curr_input);
	my $invalid_positions = "$curr_input/invalid-positions.tsv";

	my $extra_params = " --fasta $curr_input/reference.fasta ";

	if (-e "$invalid_positions")
	{
		$extra_params .= "--invalid-pos $invalid_positions";
		print "testing with invalid positions $invalid_positions\n";
	}

	my $description = `cat $curr_input/description`;
	my $expected = `cat $curr_input/expected.fasta`;
	my $expected_out_file = "$curr_input/expected.fasta";
	my $expected_positions_file = "$curr_input/expected.positions.tsv";
        my $expected_core_file = "$curr_input/expected_core.csv";

	my $vcf_dir = "$curr_input/consolidate_vcf";
  my %vcf_files = %{get_bcfs($vcf_dir)};



	test_header($curr_input);
	print "### Description ###\n";
	print "$description\n";
	print "### Expected ###\n";
	print "$expected\n";

	die "$expected_out_file does not exist" if (not -e $expected_out_file);

	my $done_testing = 0;

	my ($actual_base,$actual_out_file) = run_command(\%vcf_files,'ref',['fasta'], $extra_params);
	my $actual_positions_file = "$actual_base-positions.tsv";
        my $actual_core_file = "$actual_base-stats.csv";
	my $got = -e $actual_out_file ? `cat $actual_out_file` : 'empty file';
	print "### Got ###\n";
	print "$got\n";
	my $success = compare_files($expected_out_file,$actual_out_file);
	if ($success) #pass
	{
		pass("snv_alignment generated from data in $curr_input is valid");
	}
	$success = compare_files($expected_positions_file,$actual_positions_file);
	if ($success)
	{
		pass("positions file generated from data in $curr_input is valid");
	}

	$success = compare_files($expected_core_file,$actual_core_file);
	if ($success)
	{
		pass("core file generated from data in $curr_input is valid");
	}


	print "### done ###\n";
}

done_testing();
