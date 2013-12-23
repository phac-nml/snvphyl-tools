#!/usr/bin/env perl
# core-only.pl
# Purpose:  Given a gbk file, create a filter tab position file of all intergentic regions 
use warnings;
use strict;
use FindBin;
use lib $FindBin::Bin.'/lib';
use lib '.';
use Getopt::Long;
use Bio::SeqIO;
use autodie;

my $in_gbk;
my $out_tsv;

if (!GetOptions('i|gbk=s' => \$in_gbk,
                'out=s' => \$out_tsv)){
	die "Invalid option\n";
}

if (not $in_gbk || not -e $in_gbk){
	die "ERROR: Could not find input genbank '$in_gbk'\n";
}


if (not $out_tsv){
	die "ERROR: Was not given an output file\n";
}


open my $out,">",$out_tsv;


my $in = Bio::SeqIO->new(-format=>'genbank',-file=>$in_gbk);
my %non_coding;

while ( my $seq = $in->next_seq()){
    #determine size of genome
    my $length = $seq->length();
    my $id = $seq->display_id();
    
    my @cds;
    
    for my $feat_object ($seq->get_SeqFeatures) {
        if ( $feat_object->primary_tag eq 'CDS') {
            my ($start,$end) = ($feat_object->start,$feat_object->end);
            if ( $start > $end) {
                push @cds,{'start'=>$end,'end'=> $start};
            }
            else {
                push @cds,{'start'=>$start,'end'=> $end};
            }
        }
    }

    
    my @region;
    @cds =  sort { $a->{'start'} <=> $b->{'start'} } @cds;

    #take first one so we can simplify the loop below
    #just looking for overlapping regions
    if ( scalar @cds !=0) {
        push @region, {'start' => $cds[0]->{'start'},'end' => $cds[0]->{'end'}};
        shift @cds;
        foreach my $cds(@cds) {
            if ($region[$#region]->{'end'} > $cds->{'start'}  ) {
                $region[$#region]->{'end'} = $cds->{'end'};
            }
            else {
                push @region,{'start'=>$cds->{'start'},'end'=>$cds->{'end'}};
            }
        }
        
        
        my ($start,$end);

        if ( $region[0]->{'start'}==1) {
            $start = $region[0]->{'end'} +1;
        }
        else {
            $start = 1;
            $end = $region[0]->{'start'} -1;
            $non_coding{$id}= [([($start,$end)])];
            $start = $region[0]->{'end'} + 1;
            $end =$start;
        }
        
        shift @region;
    
        foreach my $cds (@region) {
            #record positions
            $end = $cds->{'start'}-1;
            if ( exists $non_coding{$id}) {
                push @{$non_coding{$id}},[($start,$end)];
            }
            else {
                $non_coding{$id}= [([($start,$end)])];
            }
            
            $start = $cds->{'end'} + 1
        }
        
        #ensure that we get the last possible region of the genome
        if ( $start < $length) {
            push @{$non_coding{$id}},[($start,$length)];
        }
    #check to see if the cds is the first position, if so, no range at the beginning, otherwise process it
    #no regions, the whole region is gone!
    }
    
    else {
        $non_coding{$id}= [([(1,$length)])];
    }
        


}

#print everything out
foreach my $id( keys %non_coding) {
    foreach ( @{ $non_coding{$id}} ) {
        print $out join(',' , ($id, @{$_})) . "\n";
    }
}

exit;

