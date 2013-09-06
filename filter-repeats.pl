#!/usr/bin/env perl
# filter-repeats.pl
# Purpose:  Given a nucmer repeats file, removes the "odd" repeats.

use warnings;
use strict;

my ($input) = @ARGV;

sub usage
{
	"Usage: $0 [repeats file]\n".
	"Parameters:\n".
	"\t[repeats file]:  A file containing repeats from nucmer, generated from 'show-coords -TH'\n";
}

# MAIN
die "Error: no input file defined\n".usage if (not defined $input);
die "Error: file $input does not exist" if (not -e $input);

open(my $fh, "<$input") or die "Could not open $input: $!";

while(my $line = readline($fh))
{
	chomp $line;
	my @fields = split(/\t/,$line);
	die "Error: invalid show-coords file in $input" if (@fields != 9);

	my ($start1,$end1,$start2,$end2,$length1,$length2,$pid,$contig1,$contig2) = @fields;
	my $diff_start = $start2-$start1;
	my $diff_end = $end2-$end1;

	if ($contig1 eq $contig2 and
		($diff_start == 0 or $diff_end == 0))
	{
		print STDERR "$line\n";
	}
	else
	{
		print "$line\n";
	}
}
