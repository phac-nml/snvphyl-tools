#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;
use Test::More;
use Test::Exception;
use File::Compare;
use File::Temp qw(tempdir);

my $script_dir = $FindBin::Bin;
my $rearrange_bin = "$script_dir/../rearrange_snp_matrix.pl";

my $old_env = $ENV{'PERL5LIB'};
$ENV{'PERL5LIB'} = "$script_dir/../lib:$script_dir/../cpanlib/lib/perl5:";
$ENV{'PERL5LIB'} .= $old_env if (defined $old_env);

my $command;
#==============================================================================
#=======VALID INPUT FUNCTIONALITY TESTS=============
#1 => verify that the script saves output in the proper file locations
my $create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";

system("$rearrange_bin -k -r VC-18 -s increasing -t $script_dir/tree/input/pseudoalign.phy_phyml_tree.txt -o ".$create_dir."/ -m $script_dir/tree/input/matrix.csv -p $script_dir/tree/input/pseudoalign.phy 2>&1");
ok(-e $create_dir.'/phylogeneticTree.txt', "Output tree file in correct location.");
ok(-e $create_dir.'/revisedMatrix.csv', "Revised matrix file in correct location.");

#2 => verify that the script does not make any alterations to the tree or matrix when none of the optional parameters are set
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("$rearrange_bin -t $script_dir/tree/input/pseudoalign.phy_phyml_tree.txt -o ".$create_dir."/ -m $script_dir/tree/input/matrix.csv -p $script_dir/tree/input/pseudoalign.phy 2>&1");
ok((compare($create_dir.'/phylogeneticTree.txt', $script_dir.'/tree/expected/2/phylogeneticTree.txt')==0), "Verify that the phylogenetic tree has remained unchanged.");
ok(!(-e $create_dir.'/revisedMatrix.csv'), "Verify that the revised matrix file was not generated.");

#3 => verify that the tree will convert branch lengths to total SNP estimate
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("$rearrange_bin -c -t $script_dir/tree/input/pseudoalign.phy_phyml_tree.txt -o ".$create_dir."/ -m $script_dir/tree/input/matrix.csv -p $script_dir/tree/input/pseudoalign.phy 2>&1");
ok((compare($create_dir.'/phylogeneticTree.txt', $script_dir.'/tree/expected/3/phylogeneticTree.txt')==0), "Verify that the phylogenetic tree has its branch lengths converted to total SNP estimates.");

#4 => verify that the tree is properly sorted in ascending order
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("$rearrange_bin -s increasing -t $script_dir/tree/input/pseudoalign.phy_phyml_tree.txt -o ".$create_dir."/ -m $script_dir/tree/input/matrix.csv -p $script_dir/tree/input/pseudoalign.phy 2>&1");
ok((compare($create_dir.'/phylogeneticTree.txt', $script_dir.'/tree/expected/4/phylogeneticTree.txt')==0), "Verify that the phylogenetic tree is properly sorted in increasing order.");

#5 => verify that the tree is properly sorted in descending order
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("$rearrange_bin -s decreasing -t $script_dir/tree/input/pseudoalign.phy_phyml_tree.txt -o ".$create_dir."/ -m $script_dir/tree/input/matrix.csv -p $script_dir/tree/input/pseudoalign.phy 2>&1");
ok((compare($create_dir.'/phylogeneticTree.txt', $script_dir.'/tree/expected/5/phylogeneticTree.txt')==0), "Verify that the phylogenetic tree is properly sorted in descending order.");

#6 => verify that the tree is properly re-rooted
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("$rearrange_bin -c -r VC-18 -t $script_dir/tree/input/pseudoalign.phy_phyml_tree.txt -o ".$create_dir."/ -m $script_dir/tree/input/matrix.csv -p $script_dir/tree/input/pseudoalign.phy 2>&1");
ok((compare($create_dir.'/phylogeneticTree.txt', $script_dir.'/tree/expected/6/phylogeneticTree.txt')==0), "Verify that the phylogenetic tree has its branch lengths converted to total SNP estimates.");

#7 => verify that the matrix.csv is properly re-ordered
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
system("$rearrange_bin -k -c -r VC-10 -s increasing -t $script_dir/tree/input/pseudoalign.phy_phyml_tree.txt -o ".$create_dir."/ -m $script_dir/tree/input/matrix.csv -p $script_dir/tree/input/pseudoalign.phy 2>&1");
ok((compare($create_dir.'/revisedMatrix.csv', $script_dir.'/tree/expected/4/revisedMatrix.csv')==0), "Verify that the revisedMatrix.csv is properly sorted in increasing order.");

#=========INVALID INPUT ERROR HANDLING TESTS=========
#7 => verify that the script dies properly when an invalid phylogenetic tree
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
$command = "$rearrange_bin -r VC-18 -s increasing -t $script_dir/tree/input/INVALIDpseudoalign.phy_phyml_tree.txt -o ".$create_dir."/ -m $script_dir/tree/input/matrix.csv -p $script_dir/tree/input/pseudoalign.phy 2>/dev/null";
ok(system(`$command`)!=0, "Invalid newick files are not accepted.");

#8 => verify that the script dies properly when an invalid matrix.csv file is input
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
$command = "$rearrange_bin -r VC-18 -s increasing -t $script_dir/tree/input/pseudoalign.phy_phyml_tree.txt -o ".$create_dir."/ -m $script_dir/tree/input/INVALIDmatrix.csv -p $script_dir/tree/input/pseudoalign.phy 2>/dev/null";
ok(system(`$command`)!=0, "Invalid newick files are not accepted.");

#9 => verify that the script dies properly when an invalid pseudoalign.phy file is input
$create_dir = tempdir(TEMPLATE => 'XXXXX', CLEANUP => 1) or die "Unable to create temporary file directory.";
$command = "$rearrange_bin -r VC-18 -s increasing -t $script_dir/tree/input/pseudoalign.phy_phyml_tree.txt -o ".$create_dir."/ -m $script_dir/tree/input/matrix.csv -p $script_dir/tree/input/INVALIDpseudoalign.phy 2>/dev/null";
ok(system(`$command`)!=0, "Invalid pseudoalign.phy files are not accepted.");

#=============================================================================
done_testing();
