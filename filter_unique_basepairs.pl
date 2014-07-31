#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::TreeIO;
use List::MoreUtils 'any';
use File::Basename;

main( check_inputs() );

sub main
{
    my ($tsv,$flag,$tree_file,$output, $clade_output, $snpeff_input_ref, @desired_strains) = @_;

    my @strains= map_strains($tsv, $tree_file);
    my @tsv_array = get_all_positions($tsv);


   my %snpeff_info;

    open my $outfile, '>', $output or die "Could not open file for writing: $?";

    if(scalar @desired_strains == 0)
    {
    	my %clades = make_tree_file_get_clades($tree_file, $clade_output);
    	my $index = 1;
    	foreach my $key ( keys %clades )
         {
	    my @desired;
    	    for my $strain (@{$clades{"Clade_".$index}})
    	    {
  		  push(@desired, $strain);
	    }
	    print $outfile sprintf("Clade_".$index.":,\n");
    	    process_strains(\@tsv_array, $outfile, \@desired, $flag, \@strains, $snpeff_input_ref);
    	    print $outfile sprintf("Strains:, ("."@desired".")\n\n");
    	    print $outfile "--------------------------------------------------------------------------------------------\n\n";
    	    $index++;
	 }	


        print $outfile "\n\nIndividual Genome Comparisons:,\n";
        for my $i (3.. ( scalar @strains -1 ) )
         {
            my @desired;
            push(@desired, $strains[$i]);
            print $outfile sprintf("Strains: ("."$strains[$i]".")\n");
            process_strains(\@tsv_array, $outfile, \@desired, $flag, \@strains, $snpeff_input_ref, $snpeff_input_ref);
            printf $outfile "\n\n";
            print $outfile "--------------------------------------------------------------------------------------------\n\n";
        }
    }
    else
    {
         print $outfile sprintf("Strains: ("."@desired_strains".")\n");
         process_strains(\@tsv_array, $outfile, \@desired_strains, $flag, \@strains, $snpeff_input_ref, $snpeff_input_ref);
	 printf $outfile "\n\n";
	 print $outfile "--------------------------------------------------------------------------------------------\n\n";
    }
   close $outfile;

}

sub map_strains
{
    my ($tsv, $tree_in) = @_;
    my @strains;

	open( my $pseudo_align, '<', $tsv ) or die "could not open file";

    #grab the strain names from the top of the tsv file
    my $line = <$pseudo_align>;
    chomp $line;

    @strains = split( /\s/, $line );
    my $tree = Bio::TreeIO->new(-format => 'newick', -file => $tree_in)->next_tree;

    for my $node ( $tree->get_nodes ) 
    {
        if($node->id)
        {
            my $id = $node->id;
            $id =~ s/^'//;
            $id =~ s/'$//;
    	   #it's a strain and not in strains, then we need to replace Reference with it
    	   if( ($id !~ /^-?\d+\.?\d*$/) and (index_of($id, \@strains) == -1 ) )
    	   {
               $strains[3] = $id;
    	   }
        }
	}
    return @strains;
}

sub get_all_positions
{
    my ($tsv) = @_;
    my $posn = 0;
    my @tsv_array;

    open( my $pseudo_align, '<', $tsv ) or die "could not open file";

    #grab the file headers
    my $line = <$pseudo_align>;
    chomp $line;

    while($line = <$pseudo_align>)
    {
        chomp $line;
        my @bps = split( /\s/, $line );
        push(@{$tsv_array[$posn]}, @bps);
        $posn++;
    }

    return @tsv_array;
}

sub make_tree_file_get_clades
{
    my ($tree_in, $clade_output) = @_;
	my $index = 1;

	my $tree = Bio::TreeIO->new(-format => 'newick', -file => $tree_in)->next_tree;
	my $out = new Bio::TreeIO(-file => ">$clade_output", -format => 'newick');    


    for my $node ( $tree->get_nodes ) 
    {
	
    	if($node->id)
    	{
	       if($node->id =~ /^-?\d+\.?\d*$/ )
	       {
			$node->id("Clade_".$index);
			$index++;
	       }
	       else
	       {
                        my $id = $node->id;
                        $id =~ s/^'//;
                        $id =~ s/'$//;
			#clean 's off the id
			$node->id( $id );
	       }
      }
   }

    $out->write_tree($tree);
    my %clades = get_clades_from_tree($tree);
    
    return %clades;
}

sub get_clades_from_tree
{
	my $tree = $_[0];
	my %clades;

    for my $node ( $tree->get_nodes ) 
    {
    	#clade node, strain nodes or leaves, do not have descendants
		if ( $node->id and scalar $node->get_all_Descendents > 0 )
		{
			my @strains;

			#for each node get descendants
			for my $kid ( $node->get_all_Descendents ) 
			{
				if($kid->is_Leaf)
				{
					push(@strains, $kid->id);
				}	
			}	
			push(@{$clades{$node->id}}, @strains);
		}
    }
	return %clades;
}

sub process_strains {
    my ($tsv_array, $outfile, $desired_strains_ref, $flag, $all_strains_ref, $snpeff_input_ref) = @_;
    my %strain_cols;
    my %snpeff_info;
    my @desired_strains = @{$desired_strains_ref};
    my @all_strains = @{$all_strains_ref};

    %strain_cols = find_strain_cols( \@all_strains, \@desired_strains );
    #Get any strain from this clade to be a representative
    my $strain = $desired_strains[0];
    

    print $outfile "Chromosome_Name, Position, Reference_BP, Alternate_BP";
    if ($snpeff_input_ref)
    {
         %snpeff_info = %{get_snpeff_info($snpeff_input_ref, \@desired_strains)};
         print $outfile ",Effect,Gene_Name,Old_AA/New_AA,Codon_Change,Coding,Impact,Funclass,Strand";
    }
    print $outfile "\n";

    for my $i (0.. ( scalar @$tsv_array -1 ) )
    {
        my @bps = @{@$tsv_array[$i]};
        my $valid_position;
        my $bp_match;

        #make sure the strain is valid
        if ( $flag eq "false" or $bps[2] eq "valid" ) 
        {

            #see if the given strains match at this position, Make sure that all strains in this clade all have the same bp at this location.
            $bp_match = match_strains( \@bps, \%strain_cols, $flag );
            if ($bp_match) 
            {

                #see what's going on with the other strains. Make sure no one else outside this clade has the same bp
                my $position = check_other_strains( $bp_match, \%strain_cols, @bps);
                if($position ne "")
                {
                    print $outfile $bps[0].",".$position.",".$bps[3].",".$bp_match;

                     #User provided snpeff output, so let us print the extra info
                    if ($snpeff_input_ref)
                    {
                        #Is this bp same as the reference??
                        if ($bp_match eq $bps[3])
                        {
                              print $outfile ",NO_CHANGE";
                        }
                        else
                        {
                           my $strand = "";
                           #Let us figure out the strand
                           if ($snpeff_info{$strain}{$bps[0]}{$position}{'aa'})
                           {
                              my ($old_aa, $new_aa) = split /\//, $snpeff_info{$strain}{$bps[0]}{$position}{'aa'};
                              $old_aa =~ s/[a-z]//g;
                              $new_aa =~ s/[a-z]//g;
                              
                              if (($bps[3] eq $old_aa and $bp_match eq $new_aa))
                              {
                                 $strand = "+1";
                              }
                              elsif (is_complement($bps[3], $old_aa) and is_complement($bp_match, $new_aa))
                              {
                                 $strand = "-1";
                              }
                           }

                           if (exists($snpeff_info{$strain}{$bps[0]}{$position}) and $snpeff_info{$strain}{$bps[0]}{$position}{'alt'} eq $bp_match)
                           {
                              print $outfile ",".$snpeff_info{$strain}{$bps[0]}{$position}{'effect'}.",".$snpeff_info{$strain}{$bps[0]}{$position}{'gene'}.",".$snpeff_info{$strain}{$bps[0]}{$position}{'aa'}.",".$snpeff_info{$strain}{$bps[0]}{$position}{'codon'}.",".$snpeff_info{$strain}{$bps[0]}{$position}{'coding'}.",".$snpeff_info{$strain}{$bps[0]}{$position}{'impact'}.",".$snpeff_info{$strain}{$bps[0]}{$position}{'funclass'}.",".$strand;
                           }
                           elsif ($snpeff_info{$strain}{$bps[0]}{$position}{'alt'} ne $bp_match)
                           {
                              print $outfile ",ERROR_BP_DOES_NOT_MATCH_SNPEFF_OUTPUT Ref: $snpeff_info{$strain}{$bps[0]}{$position}{'alt'} and $bp_match";
                           }
                           else
                           {
                              print $outfile ",ERROR_SNP_NOT_FOUND_IN_SNPEFF_OUTPUT";
                           }
                        }
                    }
                    print $outfile "\n";

                }
            }
         } 
    }
}

sub is_complement
{
   my ($old, $new) = @_;
   if (complement($old) eq $new)
   {
      return 1;
   }
   else
   {
      return 0;
   }
}

sub complement
{
   my $aa = shift;
   if ($aa eq "T")
   {
      return "A";
   }
   elsif ($aa eq "A")
   {
      return "T";
   }
   elsif ($aa eq "C")
   {
      return "G";
   }
   elsif ($aa eq "G")
   {
      return "C";
   }

}

sub get_snpeff_info
{
      my ($snpeff_input_ref, $desired_strains_ref) = @_;
      my %snpeff_input = %{$snpeff_input_ref};
      my @desired_strains = @{$desired_strains_ref};
      my %snpeff_info;

      my @list_of_acceptable_effects = (  "INTERGENIC", "INTERGENIC_CONSERVED", "START_GAINED", 
                                          "START_LOST", "SYNONYMOUS_START", "NON_SYNONYMOUS_START", 
                                          "NON_SYNONYMOUS_CODING", "SYNONYMOUS_CODING", "STOP_GAINED", 
                                          "SYNONYMOUS_STOP", "STOP_LOST");
      my %hash_of_acceptable_effects = map {$_ => 1} @list_of_acceptable_effects;


      foreach my $strain (@desired_strains)
      {
            next if (!(exists($snpeff_input{$strain})));
            my $file = $snpeff_input{$strain};
            
            open my $in, "<", $file or die "Could not open VCF file: ($strain: $file) $?";

            while (my $line = <$in>)
            {
                  #Skip the comments
                  next if ($line =~ /^#/);
                  chomp $line;
                  
                  #split the columns
                  my @cols = split /\t/, $line;
                  my $chromosome = $cols[0];
                  my $position = $cols[1];
                  my $reference = $cols[3];
                  my $change = $cols[4];


                  #Grab the EFF entry from the Info column in the VCF
                  my $eff_index = -1;
                  for my $i (0 .. scalar(@cols -1))
                  {
                        if ($cols[$i] =~ /EFF/)
                        {
                              $eff_index = $i;
                              last;
                        }
                  }
                  if ($eff_index == -1)
                  {
                        die "The VCF files are not properly annotated. Cannot find the EFF information. Make sure you give me the SnpEFF VCF output files\n $file: $line";
                  }

                  #Get the EFF entry
                  my @eff = split /;/, $cols[$eff_index];

                  #Split the EFF information
                  my @all_effects = split /,/, $eff[35];

                  #For now, this will only take the "acceptable effects" and there should not be multiple entries per strain.
                  #In the future, we may want to show all effects, so the information is already here. 
                  #If we have multiple entries, we'll get an "acceptable" effect to use

                  $eff_index = -1;
                  for my $x (0 .. scalar(@all_effects-1))
                  {
                        if ($all_effects[$x] =~ /EFF=/)
                        {
                              $all_effects[$x] =~ s/EFF=//;
                        }
                        my ($temp) = $all_effects[$x] =~ /(.+?)\(/;
                        if (exists($hash_of_acceptable_effects{$temp}))
                        {
                           $eff_index = $x;
                           last;
                        }
                  }
                  if ($eff_index == -1)
                  {
                        print "Error: No acceptable effect was found in the following entry in the file: $file. The first entry was arbitarily chosen for the output.\n$line\n";
                        $eff_index = 0;
                  }

                  #extract info:
                  my ($effect, $fields) = $all_effects[$eff_index] =~ /(.+?)\((.+?)\)/;
                  my @sub_fields = split /\|/, $fields;
         
                  
                  #Let's store the information we want
                  $snpeff_info{$strain}{$chromosome}{$position}{'ref'} = $reference;
                  $snpeff_info{$strain}{$chromosome}{$position}{'alt'} = $change;
                  $snpeff_info{$strain}{$chromosome}{$position}{'gene'} = $sub_fields[5];
                  $snpeff_info{$strain}{$chromosome}{$position}{'effect'} = $effect;
                  $snpeff_info{$strain}{$chromosome}{$position}{'aa'} = $sub_fields[2];
                  $snpeff_info{$strain}{$chromosome}{$position}{'codon'} = $sub_fields[3];
                  $snpeff_info{$strain}{$chromosome}{$position}{'coding'} = $sub_fields[7];
                  $snpeff_info{$strain}{$chromosome}{$position}{'impact'} = $sub_fields[0];
                  $snpeff_info{$strain}{$chromosome}{$position}{'funclass'} = $sub_fields[1];
            }

            close $file;
      }
      
      return \%snpeff_info;
}

#given an array that represents the first line of the tsv file, and an array containing
#the strains given as arguments, this sub returns a hash containg the given strain names as
#keys, and their column in the file as values
sub find_strain_cols {
    my ($all_strains, $desired_strains) = @_;
    my %strain_cols;

    foreach my $desired_strain (@$desired_strains) {
      # for (my $i = 0; $i < (scalar @$all_strains); $i++) {
        for my $i (0 .. (scalar @$all_strains -1) ) {
            my $strain = @$all_strains[$i];

            if ( $desired_strain eq $strain ) {

                #get the column position of the desired strain in the tsv file
                $strain_cols{$desired_strain} = $i;
            }
        }

        #There was no matching strain in the tsv file
        if ( !$strain_cols{$desired_strain} ) {
        	#$strain_cols{$desired_strain} = 4;	#column number from .tsv file
            die "Error: $desired_strain is not a valid strain";
        }
    }

    return %strain_cols;
}

#given an array containing the contents of a line of the tsv file, and a hash containg
#strain name/column position pairs, this sub checks the base pairs at the columns
#specified, to see if they match or not. If they do, the matching base pair letter is returned.
#If not, an empty string is returned.
sub match_strains {
    my ($bps, $strain_cols, $flag) = @_;
    my $bp_match    = "";

    foreach my $col ( values %$strain_cols ) {
    	#since we now allow for invalid rows to be included we have to make sure that we don't
    	#accidently match 2 -'s or N's so if either var eq either of those, it can't be a match
    	if(@$bps[$col] ne "-" and @$bps[$col] ne "N" and $bp_match ne "-" and $bp_match ne "N" and length(@$bps[$col]) == 1)
    	{
	        if ( $bp_match eq "") 
	        {
	            $bp_match = @$bps[$col];
	        }
	        elsif ( $bp_match ne @$bps[$col] ) {
	            $bp_match = "";
	            last;
	        }
    	}
    	else
    	{
    		$bp_match = "";
    		last;
    	}
    }

    return $bp_match;
}

#checks the other stains besides the given matching strains,
#to see if any of the other strains contain the same base pair
#as the matching ones
sub check_other_strains {
    my ($found_bp, $strain_cols, @bps) = @_;
    my @bps_ref;
    my $no_match = "true";
    my $valid_position = "";
#remove the given strains from the bps array, so we can check the remaining strains
    @bps_ref = filter_cols( $strain_cols, \@bps );

    foreach my $other_bp (@bps_ref) {
        if ( $found_bp eq $other_bp ) {
            $no_match = "";
        }
    }

    if ($no_match) {
        $valid_position = $bps_ref[1];
    }

    return $valid_position;


}

sub filter_cols {
    my ($strain_cols, $bps) = @_;
    my @col_positions;
    my $index = 0;

    #sort the column positions by reverse order,
    #so when the @bps array is spliced, the column positions
    #don't get all messed up
    foreach my $col ( values %$strain_cols ) {
        $col_positions[$index] = $col;
        $index++;
    }

    @col_positions = sort { $b <=> $a } @col_positions;

    foreach my $col (@col_positions) {
        splice( @$bps, $col, 1 );
    }

    #splice( @$bps, 3, 1 );    #splice the reference column out
    return @$bps;
}

sub index_of
{
    my ($target, $list) = @_;
	my $found = -1;
	my $index = 0;
	foreach my $item (@$list)
	{
		if($target eq $item)
		{
			$found = $index;
		}
		$index++;
	}
	return $found;
}

sub check_inputs {
    my $tsv;
    my $tree;
    my $output;
    my $flag;
    my $help;
    my $vcf_dir;
    my %vcf_input;
    my %vcf_formatted_input;
    my @strains;
    my $clade_output;

    GetOptions(
        "tsv|t=s" => \$tsv,
        "tree|p=s" => \$tree,
        "valid|v=s" => \$flag,
        "output|o=s" => \$output,
        "clade|c=s"  => \$clade_output,
        "vcf|f=s",   => \%vcf_input,
        "strains|s=s" => \@strains,
        "directory|d=s" => \$vcf_dir,
        "help|h"  => \$help
    );

    if ($help) {
        pod2usage(
            -verbose  => 99,
            -sections => "SYNOPSIS|OPTIONS|DESCRIPTION"
        );
    }

    if ( !$tsv || !$tree ) {
        pod2usage( -verbose => 1 );
    }

    if( !$output )
    {
    	$output = "output";
    }

    if ( !$clade_output)
    {
      $clade_output = "clades.tre";
    }

    if ( !$flag )
    {
    	$flag = "true";
    }

    if ($vcf_dir)
    {
      opendir(DIR, $vcf_dir) or die "Could not open VCF directory: $?";
      my @files = readdir(DIR);
      closedir(DIR);
      foreach my $file (@files)
      {
         next if($file =~ /^\./);
         my $base = basename($file);
         (my $noext = $base) =~ s/\.[^.]+$//;
         my $fullpath = $vcf_dir."/".$file;
         $vcf_formatted_input{$noext} = $fullpath;
      }
    }
    else
    {
      foreach my $key (keys %vcf_input)
      {
         my $base = basename($key);
         (my $noext = $base) =~ s/\.[^.]+$//;
         $vcf_formatted_input{$noext} = $vcf_input{$key};
      }
   }

    return $tsv, $flag, $tree, $output, $clade_output, \%vcf_formatted_input,  @strains;
}

=head1 NAME

filter_unique_basepairs.pl

=head1 VERSION

Version 1.0

=head1 SYNOPSIS

filter_unique_basepairs.pl -t tsv_file -p phylotree.tre -v <optional> -o output -s I<strain 1> I<strain 2> ... <optional>

=head1 OPTIONS

=over

=item B<-t>, B<--tsv>

The tsv file containing the pseudoalignment

=item B<-p>, B<--tree>

The .tre file that contains the data for making the tree

=item B<-v>, B<--valid>

Boolean flag regarding the validity of a row. Default is true for high confidence (do not include flag). For all matches in including non-valid (low confidence) set flag to false.

=item B<-o>, B<--output>

The file the matching positions will be written to

=item B<-f>, B<--vcf>

A VCF input file in the following format: strain_name=file_path

=item B<-d>, B<--directory>

A folder containing all of the tabular Snpeff outputs

=item B<-s>, B<--strains>

The strains you wish to find unique basepairs in

=item B<-c>, B<--clade>

The output name for the clades tree file produced

=back

=head1 DESCRIPTION

 filter_unique_basepairs determines at what positions in a pseudoalignment file for given number of strains have unique basepairs compared to the other strains in the file. 
 If you wish to compare manually selected strains be sure to enter in each strain you wish to compare

 The program will generate 2 files:
 1) a .tre files with clades labelled on the branches
 2) an output file which lists the unique base pairs found by the listed strains, if applicable, otherwise by clade.

=cut
