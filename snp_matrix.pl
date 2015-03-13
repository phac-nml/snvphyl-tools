#!/usr/bin/env perl
use strict;
use File::Basename;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Cwd;
use Getopt::Long;

my $script_name = $0;

sub usage
{
	"Usage: ".basename($script_name)." [pseudoalign.phy]\n".
	"Constructs a snp matrix from the pseudoalignment file of the pipeline\n".
	"Options:\n".
	"\t-v:  Print more information (to stderr)\n".
	"\t-o:  Prints matrix to passed file\n";
}

my $verbose = 0;
my $output = undef;
if (not GetOptions('v' => \$verbose, 'o|output=s' => \$output))
{
	die "Invalid options\n".usage;
}

my $input_file = $ARGV[0];

die "Invalid file\n".usage if (not defined $input_file);
die "File $input_file does not exist\n".usage if (not -e $input_file);
die "File $input_file is not readable\n".usage if (not -r $input_file);
die "Invalid file $input_file\n".usage if (-d $input_file);

my $out_h;

if (defined $output)
{
	open($out_h,">$output") or die "Could not open $output for writing\n".usage;
}
else
{
	$out_h = *STDOUT;
}

my $in =  new Bio::AlignIO(-file=>$input_file, -format=>"phylip",-longid=>1);
die "Could not open $input_file as phylip formatted file\n" if (not defined $in);

my $aln = $in->next_aln;
die "No phylip formatted alignment found in $input_file\n" if (not defined $aln);

my %longseq;
my @columns;
for my $seq ($aln->each_alphabetically) {
    my @nucleotides = split //, $seq->seq;
    $longseq{$seq->display_id} .= $seq->seq;
    $longseq{$seq->display_id} .= "*";
}

my @accessions = sort {$a cmp $b } keys %longseq;
for my $acc (@accessions) {
    my @chars  = split undef, $longseq{$acc};
    for (my $x=0;$x<@chars;$x++) {
	$columns[$x] .= $chars[$x];
    }   
}
my %matrix;
my $accsorter;
column: for my $column (@columns) {
    my %strain;
    my @chars = split undef, $column;
    my @accessions = sort {$a cmp $b } keys %longseq;
    for my $char (@chars) {
	my $accession = shift @accessions;
	push @{$strain{$char}}, $accession;
	
    }
    if (scalar (keys %strain)>2){ 
	warn "warning: \"$column\" has more than 2 snp differences!\n" if ($verbose); 
    }  
    if (scalar (keys %strain)<2){ 
	warn "warning: \"$column\" has no SNP differences!\n" if ($verbose); 
	next column;
    }

    my ($invalid_base) = ($column =~ /([^ACTGactg])/);
    if (defined $invalid_base)
    {
        warn "warning: \"$column\" has an invalid base ($invalid_base), skipping column\n" if ($verbose);
        next column;
    }
    
    my @nucleotides  = keys %strain;
    my %seen;
    for my $snp1 (@nucleotides) {
	snp2: for my $snp2 (@nucleotides) {
	    next snp2 if $snp1 eq $snp2; 
	    next snp2 if $seen{$snp1}{$snp2};
	    next snp2 if $seen{$snp2}{$snp1};
	    $seen{$snp1}{$snp2}++;
	    $seen{$snp2}{$snp1}++;
	    my @group1 = @{$strain{$snp1}};
	    my @group2 = @{$strain{$snp2}};
	    for my $acc1 (@group1) {
		for my $acc2 (@group2) {
		    $matrix{$acc1}{$acc2}++;
		    $matrix{$acc2}{$acc1}++;	
	}
	    }
	}
    }
}
use Data::Dumper;
my %accsorter;
my %sortseen;
# sort the accessions by values, by inverting the keys, and values
for my $accession1 (@accessions) {
	for my $accession2 (@accessions) {
		push @{$accsorter{$matrix{$accession1}{$accession2}}}, [$accession1, $accession2] unless  $sortseen{$accession2}{$accession1};
		$sortseen{$accession1}{$accession2}++;	
	}
}
my %accounter; map {$accounter{$_}++} @accessions;
my @sortedaccs;
for my $val (sort {$b <=> $a } keys %accsorter) {
	last unless keys %accounter;
	my @acclist = @{$accsorter{$val}};
	for my $accs_r (@acclist) {
		my @accs = @$accs_r;
		for my $acc (@accs) {
			if ($accounter{$acc}) {
				push @sortedaccs, $acc;
				delete $accounter{$acc};
			}	               
		}
	}
}
# print Dumper \%accsorter; exit;
print $out_h join "\t","strain",@sortedaccs, "\n";
for my $acc1 (@sortedaccs)  {
    push my @row, $acc1?$acc1:"0";
    for my $acc2 (@sortedaccs) {
	push @row, $matrix{$acc1}{$acc2}?$matrix{$acc1}{$acc2}:"0";
    }
    print $out_h join "\t", @row, "\n";
}

