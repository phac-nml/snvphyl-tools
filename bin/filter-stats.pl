#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

my ($input, $invalids, %counts, %totals, %totalsFiltered, $help);

Getopt::Long::Configure('bundling');
GetOptions(
	'i|input=s'	=> \$input,
	'a|all'		=> \$invalids,
	'h|help'	=> \$help
);
pod2usage(0) if $help;
pod2usage(1) unless $input;

#variables to track more detailed stats on variants
my $total = 0;
my $total_excluding_invalid = 0;
my $total_filtered = 0;
my $total_invalid = 0;
my $total_density = 0;
my $total_mpileup = 0;
my $total_freebayes = 0;
my $total_coverage = 0;

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
#var keeps track of the total number of N's and -'s for each chromosome 
my %total_N; 
#Go through each line and process it
while ($line = <$in>)
{
	chomp $line;
	my @entries = split(/\t/, $line);
	#Get the chromosome name which is the first thing in the line
	my $chrom = $entries[0];	

    #Increment totals
    $totals{$chrom}++;
        	
	detailed_filter_stats($entries[2]);
		
	#Not valid? Increment the total number of filtered SNP's for chromosome
	$totalsFiltered{$chrom}++ if($invalids || !($entries[2] eq "filtered-invalid"));
	#Flag to indicate whether an N or - is found in any genome for a given position:
    my $n_flag = 0;
    my $invalid_total_flag = 0;
	#Go through the genomes. First genome starts at the 4th column
	for my $i(4 .. $#entries)
	{
		my $isolate_N_flag = 1;
                if($entries[2] eq "filtered-invalid" && !$invalids){
			#we don't want to count filtered-invalid, undo the total increment done earlier
			$totals{$chrom}-- if !$invalid_total_flag;
			$invalid_total_flag = 1;
                        $isolate_N_flag = 0;
		};
		
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
		if (($entries[$i] eq "N" or $entries[$i] eq "-") && $isolate_N_flag)
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
			#Either way, we increment the total
			$counts{$chrom}{$gen}{$t}++;
			if(!$n_flag){$total_N{$chrom}++; $n_flag=1;};
		}
	}
	
}
close($in);

my ($header_out, $t_count, $t_perc, $t_totals);

foreach my $chromosome(sort {$a cmp $b} keys %counts)
{
    if(!$total_N{$chromosome}){$total_N{$chromosome} = 0;};
    #set up the header columns/rows and print the summary of combined results as the first column
	my $chromosome_total_unrounded = $totals{$chromosome} ? $total_N{$chromosome}/$totals{$chromosome} * 100 : 0;
	my $chromosome_total_percent = sprintf("%.2f", $chromosome_total_unrounded);

	$header_out = $chromosome."\t"."ALL";
	$t_count = "Total number of N's and -'s"."\t".$total_N{$chromosome};
	$t_perc = "Total percent of N's and -'s"."\t".$chromosome_total_percent;
    $t_totals = "Total number of unfiltered variants in chromosome: ".$totals{$chromosome};
	#Sort the entries by the total count of N's and -'s in descending order
	for my $genome (sort {$counts{$chromosome}{$b}{$t} <=> $counts{$chromosome}{$a}{$t} || $a cmp $b} keys %{$counts{$chromosome}})
	{	
		#Get list of genome names
		$header_out = $header_out."\t".$genome;

		#Get the counts and percentages
		my $total_count = $counts{$chromosome}{$genome}{$t};
		my $total_percent_unrounded = $totals{$chromosome} ? $total_count/$totals{$chromosome} * 100 : 0;
		my $total_percent_rounded = sprintf("%.2f", $total_percent_unrounded);

		#Concatentate this information
		$t_count = $t_count."\t".$total_count;
		$t_perc = $t_perc."\t".$total_percent_rounded;
	}
    
	#Write everything to file
	my $temp = "Chromosome\tGenomes\n".$header_out."\n".$t_count."\n".$t_perc."\n".$t_totals."\n\n";
	print $temp;

}

#print the detailed variant stats out to the screen
my $percent_filtered = ($total_filtered/$total)*100;
my $total_used = $total - $total_filtered;

if(!$invalids){
   #exclude filtered-invalid positions from output stats:	
   $total = $total - $total_invalid;
   $total_filtered = $total_filtered - $total_invalid;
   $percent_filtered = ($total_filtered/$total)*100;
   $total_used = $total - $total_filtered;	
   $total_invalid="Invalid positions not analyzed.  Please use -a flag to analyze.\n";
printf "================= Filter Summary Statistics =====================
Number of sites used to generate phylogeny: $total_used
Total number of sites identified: $total
Number of sites filtered: $total_filtered
Percentage of sites filtered: %.2f
Coverage filtered: $total_coverage
mpileup filtered: $total_mpileup
Invalid filtered: $total_invalid\n", $percent_filtered;
}
else{ 	
printf "================= Filter Summary Statistics =====================
Number of sites used to generate phylogeny: $total_used
Total number of sites identified: $total
Number of sites filtered: $total_filtered
Percentage of sites filtered: %.2f
Coverage filtered: $total_coverage
mpileup filtered: $total_mpileup
Invalid filtered: $total_invalid\n", $percent_filtered;
}

#TODO: Add the density filter stat after the vcf files are actually changed:
#Density filtered: $total_density

sub detailed_filter_stats{
    my($filter_type)= @_;
    $total++;
    if ($filter_type eq "filtered-density"){
        $total_filtered++, $total_density++;
    }
    elsif ($filter_type eq "filtered-mpileup"){
        $total_filtered++, $total_mpileup++;
    }
    elsif ($filter_type eq "filtered-coverage"){
        $total_filtered++, $total_coverage++;
    }
    elsif ($filter_type eq "filtered-invalid"){
        $total_filtered++, $total_invalid++;
    }
    elsif ($filter_type eq "filtered-freebayes"){
        $total_filtered++, $total_freebayes++;
    }
    
    return;	
}


__END__

=head1 NAME
	
	filter_stats.pl - This script prints a stat summary of the number of N's and -'s found in the psudo-alignment positions tab delimited file.

=head1 SYNOPSIS
	
	filter_stats.pl -i <input tsv file> [-a]

=head1 OPTIONS

=over 8

=item B<-i> B<--input>

The psudo-alignment positions tab delimited file

=item B<-a> B<--all>

When this option is set, the results will include all the entries marked as 'filtered-invalid'

=item B<-h> B<--help>

Prints a help message and exits

=back

=head1 DESCRIPTION

B<filter_stats.pl> prints a stat summary of the total number of N's and -'s (dashes) found in each chromosome and genome in a tab delimited file.

=cut

