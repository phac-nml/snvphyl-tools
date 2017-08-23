#!/usr/bin/env perl
use strict;
use warnings;

use FindBin;
use Bio::Phylo::IO;
use Bio::Phylo::Forest::Tree;
use Bio::TreeIO;
use Text::CSV;
use Pod::Usage;
use Getopt::Long;
use File::Temp 'tempdir';
use Pod::Usage;
use Math::Round;
use Carp qw( croak );

#=============================================================================
#Purpose:
#    Re-roots the phylogenetic tree with the indicated strain.
#
#Input: 
#    $input_taxa_tree -> Bio::Phylo::Forest::Tree phylogenetic tree to re-root
#   $newRootStrain -> String name of the strain to root with
#    $logger->reference to the Logger object.
#Output: 
#    A modified phylogenetic tree, re-rooted on the $newRootStrain
#=============================================================================
sub reRootTree
{
    my ($input_taxa_tree, $newRootStrain, $logger) = @_;
	
    my $newRoot = $newRootStrain;
    foreach my $node ($input_taxa_tree->get_nodes()){
    if($node->id && $newRoot){ 
       if($node->id eq $newRoot)
       {
          $input_taxa_tree->reroot_at_midpoint($node);
          print $logger "The phylogenetic tree has been successfully re-rooted on strain: ".$node->id()."\n";
          return;    
       }
    }
}   
    print $logger "The requested strain".$newRoot."could not be found in the phylogenetic tree.\n";
    return;
}

#==============================================================================
#Purpose:
#    Submodule to rearrange the entries in matrix.csv to match the new phylogenetic 
#ordering.
#Input: 
#    $input_taxa_tree-> Bio::Phylo::Forest::Tree phylogenetic tree to re-root
#    $inputMatrixFile-> The location of the input matrix.csv file.
#    $logger->reference to the Logger object.
#output: 
#   New matrix.csv file that matches the ordering of the re-rooted phylo tree
#==============================================================================
sub updateMatrixCsv
{
    my ($input_taxa_tree, $inputMatrixFile, $output_dir, $logger) = @_;
    
    #open a new file handle to print the output to
    open(my $revisedMatrixCsv, '>', $output_dir.'/revisedMatrix.csv') or die "ERROR: Could not open the output file: $!";
    
    #open file handle for input matrix.csv file
    open(my $data, '<', $inputMatrixFile) or die "ERROR: Could not open '$inputMatrixFile' $!\n";
    
    my $csv = Text::CSV->new({ sep_char => '\t' });
    #hash the two-dimensional matrix as 'keyColumn:keyRow' : 'value' pairs to facilitate rearrangement 
    my %matrixHash = ();
    my @strainColumn=[];
    #parse the input matrix and retain the first row for strain names in @strainNames
    while (my $input = <$data>) {
        my @line = split(/\t/x, $input);
        #add all of the strain names for each row
        if($line[0] eq 'strain'){
            @strainColumn = split(/ /, "@line");
        }
        else{
            my $strainRow = $line[0];
            my $index = 0;
            #hash the values in the matrix as 'keyColumn:keyRow' = 'value'
            foreach(@line){
                $matrixHash{$strainColumn[$index].':'.$strainRow} = $_; 
                $index++; 
            }    
        }
    }
            
    #using reference tree, print a re-ordered matrix.csv file to the indicated file handle
    print $revisedMatrixCsv "strain\t";
    foreach($input_taxa_tree->get_nodes){
        if($_->is_Leaf()){
            print $revisedMatrixCsv $_->id."\t";
        }
    }
    print $revisedMatrixCsv "\n";
    foreach( $input_taxa_tree->get_nodes) {
        if($_->is_Leaf()){
            my $rowNode = $_->id;
            print $revisedMatrixCsv $_->id."\t";
            foreach( $input_taxa_tree->get_nodes) {
                if($_->is_Leaf()){
                    my $hashQuery = $rowNode.':'.$_->id;
                    my $hashResult = $matrixHash{$hashQuery};
                    if( defined $hashResult ){ 
                    	print $revisedMatrixCsv $hashResult."\t" };
                }
            }
            print $revisedMatrixCsv "\n";          
        }
    }
    close($revisedMatrixCsv);
    close($data);
    return; 
}

#==============================================================================
#Purpose:
#    Converts the branch lengths to an estimate of the total number of SNV 
#differences 
#Input:
#    $input_taxa_tree -> Bio::Phylo::Forest::Tree phylogenetic tree to re-root
#    $inputPhyFile -> The location of the input snv_align.phy file.
#    $logger->reference to the Logger object.
#==============================================================================
sub branchLengthToSNV
{
    my ($input_taxa_tree, $inputPhyFile, $logger) = @_;
        
    #parse the total SNV's in the tree from the input .phy file:
    open(my $inputPhy, "<", $inputPhyFile) or die "ERROR: Could not open input phylogeny.phy file.";
    my $input = <$inputPhy>;
    my @line = split(/\s/x, $input);    
    
    #extract the total number of SNV's found in the third column of the first line of the phylogeny.phy file.
    my $treeTotalSNV = $line[2];
    
    my $internalNumber = 1;
    foreach my $node ( $input_taxa_tree->get_nodes ){
      my $nodeBranchLength = $node->branch_length();
      if(not defined $nodeBranchLength){
         $nodeBranchLength = 0;
         warn "Undefined node branch length set to 0.";
      }
      my $lengthToSNV = $nodeBranchLength*$treeTotalSNV;
      $lengthToSNV = sprintf "%.2f", $lengthToSNV;
      my $nodeName = $node->id();
      
      if($lengthToSNV<0.5){
          $node->branch_length(0);
      }
      else{
          #round the number printed to the closest integer value
          $node->branch_length(round($lengthToSNV));
      }
      $internalNumber++ if !$node->is_Leaf;
    }
    close($inputPhy);
    print $logger "Branch lengths have been successfully converted to total SNV's\n";
    return;
}

#==============================================================================
#Purpose:
#    Exponentially scales branch lengths on the input tree to allow the user to 
#    make the final output more easily human readable.
#Input:
#    $input_taxa_tree -> Bio::Phylo::Forest::Tree phylogenetic tree to re-root
#   $exponent -> the exponent to factor all branch lengths by
#    $logger -> Reference to a Logger object.
#Output:
#    matrix.csv file with resized branch lengths.
#==============================================================================
sub resizeTree
{
    my ($input_taxa_tree, $exponent, $logger) = @_;
    
    #**CAUTION**:depending on the exponent value, near-zero branch lengths will become unproportionately large. 
    #To avoid this, set all branch lengths that are essentially zero to a true value of zero.  This value will depend
    #on both the value of the exponent and the branch length values of the tree.
    my $minimumBranchLength = 0.009;
    #zero out the branch values that are lower than threshold to maintain tree structure for near 0 branch lengths:
    foreach( $input_taxa_tree->get_nodes){
        if($_->get_branch_length() < $minimumBranchLength){
            $_->set_branch_length(0);
            warn "Node $_->id branch length < $minimumBranchLength, setting to 0\n";
        }
    }
    
    $input_taxa_tree->exponentiate($exponent);
    print $logger "Branch lengths have been resized with an exponential factor of: ".$exponent."\n";
    return;
}

#==============================================================================
#Purpose:
#    Script to re-root a phylogenetic tree, order the tree in increasing/
#    decreasing order, convert branch lengths to total SNV's, and output a 
#    revised matrix.csv to match the phylogenetic tree.
#==============================================================================

my ($keep_tmp, $convert, $root_strain, $tree_order, $input_tree, $output_dir, $matrix_input, $input_phy, $help );

GetOptions(
    'k|keep-tmp=s' => \$keep_tmp,
    'c|convert' => \$convert,
    'r|root=s' => \$root_strain,
    's|sort=s' => \$tree_order,
    't|tree=s' => \$input_tree,
    'o|out-dir=s' => \$output_dir,
    'm|matrix=s' => \$matrix_input,
    'p|phy=s' =>    \$input_phy,
    'h|help' => \$help
);
#handle documentation and help
if ($help){
    pod2usage(    -verbose => 99,
                -sections => "SYNOPSIS|OPTIONS|DESCRIPTION");
}

pod2usage(0) unless $input_tree && $output_dir && $matrix_input && $input_phy;

#check to ensure all required command line variables are present and set default values if applicable:
$output_dir='.' if (not defined $output_dir);
die "Error: Invalid newick file." if (not -e $input_tree);
die "Error: Invalid matrix.csv file." if (not -e $matrix_input);
die "Error: Invalid snv_align.phy file." if (not -e $input_phy);

$keep_tmp = 0 if (not defined $keep_tmp);

my $job_out = tempdir('XXXXXX', CLEANUP => (not $keep_tmp), DIR => $output_dir) or croak "Could not create temp directory";
open(my $logger, '>', $job_out.'/log.txt') or croak "Could not create temporary log output.";

#parse the newick format phylogeny generated by phyml into a Bio::Phylo::Forest::Tree object
my $taxaIO = Bio::TreeIO->new(
    '-format'   => 'newick',
    '-file'   => $input_tree
);
my $taxa = $taxaIO->next_tree;

#reroot the tree if requested by user:
reRootTree($taxa, $root_strain, $logger) if defined $root_strain;
     
#parse the Bio::Tree::Tree object into a Bio::Phylo::Forest::Tree object to facilitate reordering
my $tree = Bio::Phylo::IO->parse(
       '-string' => $taxa->as_text('newick'),
       '-format' => 'newick'
)->first;

#determine whether the tree should be re-sorted in decreasing or increasing order
if(defined $tree_order){
    if($tree_order eq "decreasing"){
        $tree->ladderize();
        print $logger "The tree has been successfully sorted in decreasing order.\n";
    }
    elsif($tree_order eq "increasing"){
        $tree->ladderize(1);
        print $logger "The tree has been successfully sorted in increasing order.\n";
    }
}
#create a matrix.csv file to reflect the changes made to the phylogenetic tree
if ((defined $root_strain) || (defined $tree_order)) { updateMatrixCsv($tree, $matrix_input, $output_dir, $logger) };
#convert branch lengths to total SNV estimate
if (defined $convert){ branchLengthToSNV($tree, $input_phy, $logger) };
     
#print the final newick formatted tree to a file that can be opened by any tree viewing program
open(my $treeout, '>', $output_dir.'/phylogeneticTree.txt') or croak "Could not open output file: $!";

#SMELLY: hack to remove any quotes that are added (for unknown reasons) to the node labels by the Bio::Phylo::IO module
my $quotes_removed = $tree->to_newick();
$quotes_removed =~ s/'//xg;

#print $parentheses_removed; 
print $treeout $quotes_removed;
close($treeout);
exit;
=head1 NAME

rearrange_snv_matrix.pl - Script to re-root a phylogenetic tree, order the tree in increasing/decreasing order, convert branch lengths to total SNV's, and output a revised matrix.csv to match the phylogenetic tree. 

=head1 SYNOPSIS

 rearrange_snv_matrix.pl -t input_tree -o output_dir -m matrix.csv -p snv_align.phy

=head1 OPTIONS

=over

=item B<-t>, B<--tree> [required]

Newick input file describing the phylogenetic tree.

=item B<-o>, B<--out-dir> [required]

The directory for output files.

=item B<-p>, B<--phy> [required]

Input snv_align.phy file. 

=item B<-m>, B<--matrix> [required]

Input matrix.csv file.

=item B<-k>, B<--keep-tmp> [optional]

Keep the temp log upon exiting the script.

=item B<-c>, B<--convert> [optional]

Convert the branch lengths to an estimate of the total SNV number.

=item B<-r>, B<--root> [optional] 

The name of the strain to use as the root for the phylogenetic tree.

=item B<-s>, B<--sort> [optional]

Either 'increasing' or 'decreasing', indicating the order in which to sort nodes in the output phylogenetic tree.

=item B<-h>, B<--help> [optional]

To display help message

=back

=head1 DESCRIPTION

rearrange_snv_matrix steps:

1.  Takes a newick formatted phylogenetic tree and reroots the tree on the strain indicated by --root.  The tree can then be sorted in increasing or decreasing order and branch lengths are converted and displayed as the total number of SNV's.  The new phylogenetic tree is output in newick format as a text file named phylogeneticTree.txt.

2.  A new matrix.csv file is generated that matches the ordering of the phylogenetic tree done in step 1 above.

=cut
