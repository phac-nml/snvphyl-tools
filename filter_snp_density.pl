#!usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use Pod::Usage;


__PACKAGE__->run unless caller;

sub usage{
   print "Invalid options in script parameters. Please use format:";
   print "filter_snp_density.pl --bcf1=[PATH/TO/bcf1.bcf] --bcf2=[PATH/TO/bcf2.bcf] --bcfX=[PATH/TO/BCFX.bcf] --threshold [density threshold value] --output-file [PATH/TO/OUTPUT.txt]"
}

sub run{
my (%bcf_files, $output_file, $cores, $density_threshold, $region_output, $help);

if (!GetOptions('bcf=s' => \%bcf_files,
                    'output-file=s' => \$output_file,
                    'region-file=s' => \$region_output,
                    'cores=i' => \$cores, 
                    'threshold=s' => \$density_threshold,
                    'help|h' => \$help)){
        die "Invalid option\n".usage;
    }

   #Perform the input validation here:

   my @files = values %bcf_files;
   if(not defined $cores){$cores=1};
   if (@files <= 0){ die "Error: No bcf files input." };

   foreach(@files){
      if(!-e $_){ die "Error: Unable to find bcf file."};
   }

   die "No output location provided." if (not defined $output_file);
   
   my @results = analyze_bcf_density(\%bcf_files, $density_threshold, $cores, $output_file, $region_output);

   print @results;
}

sub analyze_bcf_density{

   my ($files, $density_threshold, $cores, $output_file, $region_output) = @_;
   my $pm = new Parallel::ForkManager($cores);
   my @results;
  
   $pm->run_on_finish(
      sub{
         my($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
         if(defined($data_structure_reference)){
            push @results, ${$data_structure_reference};
         }
         else{
            die "Did not receive data from one of the children."
         }
      }
   );

   foreach my $file (keys %$files){
      $pm->start && next;
      my $input_file = $files->{$file};
      
      my $result = `bcftools plugin filter_snp_density $input_file -O b -o $output_file b -- --filename $input_file --region_file $region_output --threshold $density_threshold`;
      die "Error: filter SNP density plugin exited with error while working with $input_file.\n" if (not defined $result); 
     
      $pm->finish(0, \$result);
   }
   $pm->wait_all_children;
   return @results; 
}

=pod

=head1 NAME

filter_snp_density.pl - Script for filtering high density variants within the bcf files provided.  Outputs a tab seperated list describing the high density regions within a genome.

=head1 VERSION

This documentation refers to filter_snp_density.pl version 0.0.1

=head1 SYNOPSIS

filter_snp_density.pl --bcf1=[PATH/TO/bcf1.bcf] --bcf2=[PATH/TO/bcf2.bcf] --bcfX=[PATH/TO/BCFX.bcf] --threshold [density threshold value] --output-file [PATH/TO/OUTPUT.txt] --output-regions [PATH/TO/OUTPUT_REGIONS.tsv]

=head1 OPTIONS

=over

=item B  
