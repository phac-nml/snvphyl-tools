#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Time::localtime;

print_to_file(get_inputs());

sub print_to_file
{
	my $tsv = $_[0];
	my $vcf_dir = $_[1];
	my $dest = $_[2];

	#grab the headers from the tsv to get the genome names
	open(my $pseudo_align, '<', $tsv) or die ("Cannot open file $tsv");
	my $headers = <$pseudo_align>;
	my @header_cols = split(/\s/, $headers);

	my %all_bps;

	#remove unneeded cols
	for (0 .. 3)
	{
		shift @header_cols;
	}

	#go through tsv file and build hash of hash table
	while (my $line = <$pseudo_align>)
	{
		chomp $line;
		my @base_pairs = split (/\s/, $line);

		my $position = $base_pairs[1];
		my $status = $base_pairs[2];

		#only take valid reads
		if ($status eq "valid")
		{
			#don't need the first few cols
			for (0 .. 3)
			{
				shift @base_pairs;
			}

			#for each genome, put the base pair for the current position in the appropriate hash
			for (my $i = 0; $i < scalar @header_cols; $i++)
			{
				my $genome = $header_cols[$i];
			
				$all_bps{$genome}{$position} = $base_pairs[$i];
			}
		}
	}
	close ($pseudo_align);

	#go through each genome's vcf file, look for position matches and steal the base pair sequences for the new vcf
	foreach my $genome (keys %all_bps)
	{
		open(my $new_vcf, '>', "$dest\/".$genome.".vcf") or die ("Cannot open file $dest\/$genome.vcf"); 

                print $new_vcf "##fileformat=VCFv4.1\n";
                print $new_vcf "##fileDate=".get_date()."\n";
                print $new_vcf "##commandline=\"".(join " ", $0, "-t", $tsv)."\"\n";
                print $new_vcf "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">\n";
                print $new_vcf "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
                print $new_vcf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n";
	
		open(my $freebayes_file, "<", $vcf_dir."/".$genome.".vcf") or die ("Cannot open file $genome.vcf");
	
		#go through each position in the fb file
		while (my $fb_line = <$freebayes_file>)
		{
			chomp $fb_line;
			my @fb_cols = split(/\s+/, $fb_line);
		
			#find where the actual data in the fb file is	
			next unless ($fb_cols[1] && $fb_cols[1] =~ /\d/);

			#see if the positions in the freebayes file are in the hash
			if ($all_bps{$genome}{$fb_cols[1]})
			{
				print $new_vcf $fb_cols[0]."\t".$fb_cols[1]."\t.\t".$fb_cols[3]."\t".$fb_cols[4]."\t50000\t.\tDP=1000\n";
			}
		}
		close($freebayes_file);
		close($new_vcf);
	}
}

sub get_date
{
	my $tm = localtime;
	my $year = 1900+$tm->year;
	my $month = "0".($tm->mon + 1);
	my $formatted_date = $year.$month.$tm->mday;
	return $formatted_date;
}

sub get_inputs
{
	my $tsv;
	my $vcf_dir;
	my $dest;

	GetOptions ("tsv|t=s" => \$tsv,
	    "vcf_dir|v=s" => \$vcf_dir,
	    "dest|d=s" => \$dest);

	if (!$tsv || !$vcf_dir || !$dest)
	{
		pod2usage(-verbose => 0);
	}
	
	if (!(-e $tsv))
	{
		die ("$tsv is not a file");
	}
	if (!(-d $vcf_dir))
	{
		die ("$vcf_dir is not a directory");
	}
	if (!(-r $tsv) || !(-r $vcf_dir))
	{
		die ("You do not have proper permissions");
	}
	if (!(-e $dest))
	{
		mkdir $dest;
	}

	return ($tsv, $vcf_dir, $dest);
}

=head1 NAME

tsvToVcf.pl

=head1 VERSION

Version 1.0

=head1 SYNOPSIS

tsvToVcf.pl -t pseaudoalignment_tsv_file -v freebayes_vcf_directory -d destination directory

Options:
-t
	The tsv file containing the pseudoalignment

-v
	The directory containing vcf files from freebayes

-d
	The directory the new vcf files will be output to

=cut

