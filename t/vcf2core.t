#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;
use File::Temp 'tempdir';
use Getopt::Long;

my $script_dir = $FindBin::Bin;

my $old_env = $ENV{'PERL5LIB'};
$ENV{'PERL5LIB'} = "$script_dir/../../../lib:$script_dir/../../../cpanlib/lib/perl5:";
$ENV{'PERL5LIB'} .= $old_env if (defined $old_env);

my $vcf_align_bin = "$script_dir/../vcf2core.pl";

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
    my ($pileup_dir,$reference,$coverage_cutoff,$positions,$extra_params,$pileup_vcfs) = @_;

	my ($actual_out_base) = tempdir('CLEANUP'=>1);
	$extra_params = '' if (not defined $extra_params);
	my @out_files = ();
       push @out_files,"$actual_out_base/results.csv";
    
    
	my $command;
        if ($pileup_vcfs) {
            my $singles .= " --mpileup " . join (" --mpileup " , map { " $_=" . $pileup_vcfs->{$_} } keys %$pileup_vcfs);
            $command = "$vcf_align_bin $singles --fasta $reference  --coverage-cutoff $coverage_cutoff --output-base $actual_out_base  --positions $positions  $extra_params > $out_files[0] 2> /dev/null";
        }
        else {
            $command = "$vcf_align_bin --mpileup-dir $pileup_dir --fasta $reference --output-base $actual_out_base --positions $positions --coverage-cutoff $coverage_cutoff $extra_params > $out_files[0] 2> /dev/null";
        }

	print "## Running $command\n\n";
	system($command) == 0 or die "Could not run command $command: $!";

	return ($actual_out_base,@out_files);
}

sub get_vcfs
{
    my ($in_dir) = @_;
    my %vcf;
    my $dh;
    
    opendir($dh, $in_dir) or die "error opening directory $in_dir: $!";
    %vcf = map { /^(.*)\.vcf\.gz$/; $1 => "$in_dir/$_"} grep { /\.vcf\.gz$/ } readdir($dh);
    closedir($dh);

    return \%vcf
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
my @in_files = sort grep {$_ !~ /^\./} readdir($in_h);
closedir($in_h);

print "Testing all input variants in $input_dir\n";
for my $dir (@in_files)
{
	my $curr_input = "$input_dir/$dir";
	next if (not -d $curr_input);
	my $invalid_positions = "$curr_input/invalid-positions.tsv";

	my $extra_params = '';
	if ($dir !~ /^noN/)
	{
		$extra_params .= '';
	}

	if (-e "$invalid_positions")
	{
		$extra_params .= "--invalid_pos $invalid_positions";
		print "testing with invalid positions $invalid_positions\n";
	}

	my $description = `cat $curr_input/description`;
	my $positions_file = "$curr_input/expected.positions.tsv";
        my $expected_out_file = "$curr_input/expected_core.csv";
        my $ref = "$curr_input/reference.fasta";
	my $pileup_dir = "$curr_input/pileup";
        my %pileup_vcfs = %{get_vcfs($pileup_dir)};
        

        
	test_header($curr_input);
#	print "### Description ###\n";
#	print "$description\n";

	die "$expected_out_file does not exist" if (not -e $expected_out_file);

	my $done_testing = 0;

	my ($actual_base,$actual_out_file) = run_command($pileup_dir,$ref,$coverage_cutoff,$positions_file, $extra_params);
	my $actual_positions_file = "$actual_base-positions.tsv";
	my $got = `cat $actual_out_file`;
	print "### Got ###\n";
	print "$got\n";
	my $success = compare_files($expected_out_file,$actual_out_file);
	if ($success) #pass
	{
		pass("pseudoalignment generated from data in $curr_input is valid");
	}


	test_header($curr_input);
#	print "### Description ###\n";
#	print "$description\n";

	die "$expected_out_file does not exist" if (not -e $expected_out_file);


        ($actual_base,$actual_out_file) = run_command($pileup_dir,$ref,$coverage_cutoff,$positions_file, $extra_params,\%pileup_vcfs);
        $actual_positions_file = "$actual_base-positions.tsv";
        $got = `cat $actual_out_file`;
	print "### Got ###\n";
	print "$got\n";
	$success = compare_files($expected_out_file,$actual_out_file);
	if ($success) #pass
	{
		pass("pseudoalignment generated from data in $curr_input is valid");
	}
        
	print "### done ###\n";
        
        
}

done_testing();
