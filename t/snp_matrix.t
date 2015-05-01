#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;
use Test::More;

my $script_dir = $FindBin::Bin;

my $old_env = $ENV{'PERL5LIB'};
$ENV{'PERL5LIB'} = "$script_dir/../lib:$script_dir/../cpanlib/lib/perl5:";
$ENV{'PERL5LIB'} .= $old_env if (defined $old_env);

my $matrix_bin = "$script_dir/../snp_matrix.pl";

my $matrix_dir = "$script_dir/matrix";
my $input_dir = "$matrix_dir/input";
my $good_out_dir = "$matrix_dir/output";

opendir(my $in_h,$input_dir) or die "Could not open $input_dir: $!";
my @in_files = sort {$a =~ /-(\d+)\.phy/; my $x = $1; $b =~ /-(\d+)\.phy/; my $y = $1; $x <=> $y} grep {$_ !~ /^\./} readdir($in_h);
closedir($in_h);

print "Testing all input pseudoalignments in $input_dir\n";
for my $file (@in_files)
{
	print "\n### Testing $input_dir/$file ###\n";

	my $done_testing = 0;
	my $good_out_file = "$good_out_dir/$file.out";
	my $command = "$matrix_bin $input_dir/$file 2>/dev/null";
	
	my $actual_result = `$command`;
	my @actual_lines = split(/\n/,$actual_result);

	open(my $out_h, $good_out_file) or die "Could not open $good_out_file: $!";
	while(not $done_testing and (defined (my $expected_line = readline($out_h))))
	{
		chomp $expected_line;
		my $actual_line = shift(@actual_lines);
		if ($actual_line ne $expected_line)
		{
			is($actual_line,$expected_line,"lines \"$actual_line\" and \"$expected_line\" differ");
			$done_testing = 1;
		}
	}
	if (not $done_testing) #pass
	{
		pass("SNP Matrix generated from $input_dir/$file is valid");
	}
	print "### done ###\n";
	close($out_h);
}

done_testing();
