package Reporter;

use strict;
use warnings;

sub new{
    my ($class) = @_;
    
    my $self = {};
    bless($self, $class);
    
    return $self;
}

#------------------------------------------------------------------------
#Purpose:
#    Analyzes the results from the verify_mapping.pl script to determine 
#the quality with which reads have mapped to the reference strain.  Will
#mark the analysis as a PASS, FAIL, or WARN and log additional info
#if relevant.
#Output:
#   Formatted text output for the log file, and json data.
#------------------------------------------------------------------------
sub record_read_mapping{
	
    my($self, @bams) = @_;
    #record thresholds used for depth of coverage (default = 15x)
    #parse through the % mapped for each ba file:
    my $bam_string;
	my $x = 1;
    foreach(@bams){
	    $bam_string .= '--bam bam'.$x.'=$_ ';
	}
	my $mapping_results = `perl ../verify_mapping_quality.pl $bam_string --min-map 90 --min-depth 15`;
	
	return $mapping_results;
	#if any of the reads map to less than 80% of the reference, FAIL the run
	
	#if any of the reads map >80% but <90%, add a WARNING flag
	
	#otherwise, everything is good, PASS the run 
	
}

#------------------------------------------------------------------------
#Purpose:
#    Display the output from the filter-stats.pl script to the user in an
#easily readable format for the reporter log file.  Converts the data into 
#json format as well.
#------------------------------------------------------------------------
sub record_filter_stats{
	#display total number of SNP's used to generate the phylogeny
	#display total number of core SNP's found before filtering
	#display total number and percentage of SNP's filtered
	#display the number and percentage of SNP's filtered, for each filter type
	
	#FUTURE TASK: based on organism type (S. enteriditis, Campylobacter, Ecoli, etc)
	#generate a warning if the number of SNP's used to generate the 
	#phylogenies is too low
	
}

#------------------------------------------------------------------------
#Purpose:
#   Prints user-provided information about the reference.  Converts the 
#information into json format.
#------------------------------------------------------------------------
sub record_reference_info{
    #identifier
    #type of sequencer used to generate the reference
    #are plasmids present/position of plasmids
    #is denovo assembly?
    	#YES: N50, number of contigs, coverage, min/max contig lengths		
}

#------------------------------------------------------------------------
#Purpose:
#    Record the size of all read/alignment files found in the pipeline/
#------------------------------------------------------------------------
sub record_file_sizes{
    #record all file sizes for each bam file and report any small files
    #as a quality fail
    #print "File\tSize";
    #foreach(%bam_keys){
	#	my $result = `ls -l $_`;
	#	print "$_\t$result";
	#}
	#use the average file size to determine if any of the files are far
	#smaller than they should be.  Alternatively, identify very large
	#input files that may need to be downsized for the run.
		
}

#------------------------------------------------------------------------
#Purpose:
#    Record all of the parameters used for the various steps in the 
#pipeline.  All information required to replicate the pipeline run.
#------------------------------------------------------------------------
sub record_run_parameters{
	#Files used
	#All parameter values
}

#------------------------------------------------------------------------
#Purpose:
#    Reports the gene within which a particular variant occurs and the 
#type of change that the variant will cause to the nucleotide sequence
#of the resultant protein, if applicable.
#------------------------------------------------------------------------
sub record_snp_eff{
		
}

1;

