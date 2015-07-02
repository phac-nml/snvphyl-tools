package Reporter;

use strict;
use warnings;

use String::Util 'trim';
use JSON;

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
sub record_read_mapping_data{

    my($self, @bams) = @_;
    
    #the variables to set that will be recorded
    my $reference_size;
    my @warnings;
    my $status = 'PASSED';
    my $common_problems;
    my $common_solutions;
    my $output_data={};
    my $depth;
    
    #record text of why poor mapping can cause pipeline problems
    $common_problems = "Poorly mapped reads will greatly reduce the number of\n 
                        usable SNP's and produce a final tree of lesser quality.\n 
                        Poor mapping occurs for several reasons, including:\n
                        1. Insufficient amount of sequencing data\n
                        2. Sample contamination.\n
                        3. Imporperly identified/labelled samples.\n
                        4. Poor reference strain used.\n";    
    
    $common_solutions = "Depending on the cause for poor mapping, several potential solutions exist:\n
                         1. Top up the sequencing data for the poorly mapped isolates.\n
                         2. Remove isolates that are mapping poorly from the analysis.\n
                         3. Ensure the samples are identified and labelled correctly.\n
                         4. Select a different reference strain to map to.\n";
                         
    #record thresholds used for depth of coverage (default = 15x)
    #parse through the % mapped for each ba file:
    my $bam_string;
	my $x = 1;
    foreach(@bams){
	    $bam_string .= '--bam bam'.$x.'=$_ ';
	}
	my $mapping_results = `perl ../verify_mapping_quality.pl $bam_string --min-map 90 --min-depth 15`;
	
	#parse the results from the verify_mapping_quality.pl script
	for(split /^/, $mapping_results){	
		if($_ =~ 'MINIMUM DEPTH'){
		  	$depth = trim((split ':', $_)[1]);
		}
		if($_ =~ 'REFERENCE GENOME'){
			$reference_size = trim((split ':', $_)[1]);
		}
		if($_ =~ '%'){
			push @warnings, $_;
		}
	}

	#determine if the QC check should FAIL or throw a WARNING
	foreach(@warnings){
		my $percent_mapped = trim((split ':', $_)[1]);
		$percent_mapped = trim((split '%', $percent_mapped)[0]);
		if($percent_mapped < 80){
			$status = 'FAILED';
		}
		elsif(($percent_mapped <= 90) && ($status eq 'PASSED')){
			$status = 'WARNING';
		}
	}
	
	#set the output values for the map
    $output_data->{'size'} = $reference_size;
    $output_data->{'depth'} = $depth;
    $output_data->{'status'} = $status;
    $output_data->{'Problem strains'} = \@warnings;
    $output_data->{'Common problems'} = $common_problems;
    $output_data->{'Common solutions'} = $common_solutions;
    
	return $output_data;
}

#------------------------------------------------------------------------
#Purpose:
#    Display the output from the filter-stats.pl script to the user in an
#easily readable format for the reporter log file.  Converts the data into 
#json format as well.
#------------------------------------------------------------------------
sub record_filter_stats{
	
	my($self, $pseudoalign_fp) = @_;
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
	
	my ($reference_id, $sequencer_type, $source, $assembly, $size);
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
	
	my($self, @files) = @_;
	
	my @fileSizes;
	my $min_file_size;
	my $max_file_size;
	
	foreach(@files){
	    my $output = `ls -l --block-size=M $_`;
	    my $size = (split " ", $output)[4];
	    push @fileSizes, $size;	
	}
	
	
	
	foreach(@fileSizes){
		if($_ < $min_file_size){
			#flag the file
		}
		elsif($_ > $max_file_size){
			#flag the file as requiring downsampling
		}
	}
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
	
	my($self) = @_;
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

