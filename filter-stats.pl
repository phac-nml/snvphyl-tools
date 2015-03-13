#!/usr/bin/env perl
#filer-stats.pl

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

my ($input, $output, $invalids, %counts, %totals, $help);

Getopt::Long::Configure('bundling');
GetOptions(
	'i|input=s'	=> \$input,
	'o|output=s'	=> \$output,
	'a|all'		=> \$invalids,
	'h|help'	=> \$help
);
pod2usage(1) if $help;
pod2usage(1) unless $input && $output;


#These variables are just used for the hash to store different information
#n stores the number of N's found in that genome
my $n = 0;
#d stores the number of -'s (dashes) found in that genome. Can be implemented later
my $d = 1;
#t stores the total number of N's and -'s. Can be implemented later
my $t = 2;


#Open file
open my $in, "<", $input or die "Could not open $input!";

my $line = <$in>;
chomp $line;
#Get the list of genome names from the header. 
my @header = split (/\t/, $line);

#Go through each line and process it
while ($line = <$in>)
{
	chomp $line;
	my @entries = split(/\t/, $line);
	#Get the chromosome name which is the first thing in the line
	my $chrom = $entries[0];

	#Does the user want the filtered-invalid entries included? If not, skip
	next if ($entries[2] eq "filtered-invalid" and !$invalids);

        #Increment totals
        $totals{$chrom}++;

	#Valid? No point in doing all the work. Skip
	next if ($entries[2] eq "valid");


	#Go through the genomes. First genome starts at the 4th column
	for my $i(4 .. $#entries)
	{
		#Get the name of the current genome we're looking at
		my $gen = $header[$i];

		#Make the hash entry if it doesn't exit
		if (!(exists($counts{$chrom}{$gen}{$n})))
		{
			$counts{$chrom}{$gen}{$n} = 0;
		}
		if (!(exists($counts{$chrom}{$gen}{$d})))
		{
			$counts{$chrom}{$gen}{$d} = 0;
		}
		if (!(exists($counts{$chrom}{$gen}{$t})))
		{
			$counts{$chrom}{$gen}{$t} = 0;
		}
		

		#If entry is N or -, then go count things
		if ($entries[$i] eq "N" or $entries[$i] eq "-")
		{
			#N
			if ($entries[$i] eq "N")
			{
				$counts{$chrom}{$gen}{$n}++;
			}
			#-
			elsif ($entries[$i] eq "-")
			{
				$counts{$chrom}{$gen}{$d}++;
			}

			#Either way, we increment he total
			$counts{$chrom}{$gen}{$t}++;
		}
	}
	
}
close($in);


#Now let's start writing the output file.
open my $out, ">", $output or die "Could not open $output for writing!";

my ($header_out, $t_count, $t_perc);

foreach my $chromosome(sort {$a cmp $b} keys %counts)
{
	#set up the header columns/rows
	$header_out = $chromosome;
	$t_count = "Total number of N's and -'s";
	$t_perc = "Total percent of N's and -'s";
		
	#Sort the entries by the total count of N's and -'s in descending order
	for my $genome (sort {$counts{$chromosome}{$b}{$t} <=> $counts{$chromosome}{$a}{$t}} keys %{$counts{$chromosome}})
	{	
		#Get list of genome names
		$header_out = $header_out."\t".$genome;

		#Get the counts and percentages
		my $total_count = $counts{$chromosome}{$genome}{$t};
		my $total_percent_unrounded = $total_count/$totals{$chromosome} * 100;
		my $total_percent_rounded = sprintf("%.2f", $total_percent_unrounded);

		#Concatentate this information
		$t_count = $t_count."\t".$total_count;
		$t_perc = $t_perc."\t".$total_percent_rounded;
		
	}

	#Write everthing to file
	my $temp = "Chromosome\tGenomes\n".$header_out."\n".$t_count."\n".$t_perc."\n\n";
	print $out $temp;

}

close($out);




__END__

=head1 NAME
	
	filter_stats.pl - This script prints a stat summary of the number of N's and -'s found in the psudo-alignment positions tab delimited file.

=head1 SYNOPSIS
	
	filter_stats.pl -i <input tsv file> -o <output file name> [-a]

=head1 OPTIONS

=over 8

=item B<-i> B<--input>

The psudo-alignment positions tab delimited file

=item B<-o> B<--output>

The name of the output file

=item B<-a> B<-all>

When this option is set, the summary will include all the entries marked as 'filtered-invalid'

=item B<-h> B<--help>

Prints a help message and exits

=back

=head1 DESCRIPTION

B<filter_stats.pl> prints a stat summary of the total number of N's and -'s (dashes) found in each chromosome and genome in a tab delimited file.

=cut

