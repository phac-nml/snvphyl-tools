package CompareFiles;
# Used to compare two different files for testing

use File::Copy;
use File::Basename;
use Test::Deep;

use strict;
use warnings;

sub compare_phylip_files
{
	my ($file1, $file2) = @_;

	die "error: file1 not defined" if (not defined $file1);
	die "error: file2 not defined" if (not defined $file2);
	die "error: $file1 does not exist" if (not -e $file1);
	die "error: $file2 does not exist" if (not -e $file2);

	my $data1 = _read_phylip_file($file1);
	my $data2 = _read_phylip_file($file2);

	return eq_deeply($data1, $data2);
}

sub _read_phylip_file
{
	my ($file) = @_;

	my %data;
	open(my $fh, "<$file") or die "Could not open $file";
	my $header_line = readline($fh);
	chomp $header_line;
	$data{'header'} = $header_line;

	while(my $line = readline($fh))
	{
		chomp $line;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		next if ($line eq "");
		push(@{$data{'lines'}}, $line);
	}

	close($fh);

	return \%data;
}

1;
