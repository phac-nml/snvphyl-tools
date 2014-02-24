#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

scan_file(check_inputs());

sub scan_file
{
	my $tsv = $_[0];
	my @desired_strains = @ARGV;
	my @all_strains;
	my %strain_cols;

	open (my $pseudo_align, '<', $tsv);

	#grab the strain names from the top of the tsv file
	my $line = <$pseudo_align>;
	chomp $line;
	@all_strains = split(/\s/, $line);

	#find the column positions of the strain names given as arguments
	%strain_cols = find_strain_cols(\@all_strains, \@desired_strains);

	#go through each line of the file, looking for positions where
	#the given strains differ from ervery other strain
	while ($line = <$pseudo_align>)
	{
		chomp $line;

		my $got_a_match = "true";
		my $valid_position;
		my @bps = split(/\s/, $line);
		my $bp_match;
	
		#make sure the strain is valid
		if ($bps[2] eq "valid")
		{
			#see if the given strains match at this position
			$bp_match = match_strains(\@bps, \%strain_cols);
					
			if ($bp_match)
			{
				#see what's going on with the other strains
				check_other_strains($bp_match, \%strain_cols, @bps);
			}
		}
	}
}

#given an array that represents the first line of the tsv file, and an array containing
#the strains given as arguments, this sub returns a hash containg the given strain names as
#keys, and their column in the file as values
sub find_strain_cols
{
	my $all_strains = $_[0];	
	my $desired_strains = $_[1];
	my %strain_cols;	

	foreach my $desired_strain (@$desired_strains)
	{
		for (my $i = 0; $i < scalar @$all_strains; $i++)
		{
			my $strain = @$all_strains[$i];
			
			if ($desired_strain eq $strain)
			{
				#get the column position of the desired strain in the tsv file
                                $strain_cols{$desired_strain} = $i;
			}
		}

		#There was no matching strain in the tsv file
		if (!$strain_cols{$desired_strain})
		{
			die "Error: $desired_strain is not a valid strain"
		}
	}

	return %strain_cols;
}

sub match_strains
{
	my $bps = $_[0];
	my $strain_cols = $_[1];
	my $bp_match = "";

	foreach my $col (values %$strain_cols)
        {
	        if ($bp_match eq "")
                {
	                $bp_match = @$bps[$col];
                }
                elsif ($bp_match ne @$bps[$col])
                {
			#think something buggy is going on here
	                $bp_match = "";
			last;
                }
        }

	return $bp_match;
}

sub check_other_strains
{
	my $found_bp = $_[0];
	shift @_;
	my $strain_cols = $_[0];
	shift @_;
	my @bps = @_;
	my $bps_ref;
	my $no_match = "true";

	$bps_ref = filter_cols($strain_cols, \@bps);

	foreach my $other_bp(@$bps_ref)
	{
		if ($found_bp eq $other_bp)
		{
			$no_match = "";
		}
	}

	if ($no_match)
	{
		my $valid_position = @$bps_ref[1];
                print $valid_position."\n";
        }

}

sub filter_cols
{
	my $strain_cols = $_[0];
	my $bps = $_[1];
	my @col_positions;
	my $index = 0;

	#sort the column positions by reverse order,
	#so when the @bps array is spliced, the column positions
	#don't get all messed up
	foreach my $col (values %$strain_cols)
	{
		$col_positions[$index] = $col;
		$index++;
	}

	@col_positions = sort {$b <=> $a} @col_positions;

	foreach my $col (@col_positions)
	{
		splice (@$bps, $col, 1);
	}

	splice (@$bps, 3, 1); #splice the reference column out

	return $bps;
}

sub check_inputs
{
	my $tsv;

	GetOptions ("tsv|t=s" => \$tsv);

	if (!$tsv)
	{
		pod2usage(-verbose => 0);
	}

	if (scalar @ARGV < 1)
	{
		die "Error: At least one strain must be provided as an argument";
	}		

	return $tsv;
}

=head1 NAME

filter_unique_basepairs.pl

=head1 VERSION

Version 1.0

=head1 SYNOPSIS

filter_unique_basepairs.pl -t tsv_file I<strain 1> I<strain 2> ...

-t

        The tsv file containing the pseudoalignment

strains

        The strains you wish to find unique basepairs in

=cut
