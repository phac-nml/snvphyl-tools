package Reporter;

use strict;
use warnings;

use String::Util 'trim';
use JSON;
use Switch;

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
    #parse through the % mapped for each bam file:
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
	
	my $result = `perl ../filter-stats.pl -i $pseudoalign_fp`;
	my $filter_stats = {};
	
    for(split /^/, $result){
    	switch($_){
    		case {$_ =~ 'generate phylogeny'} {$filter_stats->{'total_used'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'sites identified'} {$filter_stats->{'total_unfiltered'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'sites filtered'} {$filter_stats->{'sites_filtered'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'Coverage filtered'} {$filter_stats->{'filtered_coverage'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'mpileup filtered'} {$filter_stats->{'filtered_mpileup'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'Density filtered'} {$filter_stats->{'filtered_density'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'Invalid filtered'} {$filter_stats->{'filtered_invalid'}=trim((split ":", $_ )[1])}
        }
    }    
    print to_json($filter_stats);
}

#------------------------------------------------------------------------
#Purpose:
#   Prints user-provided information about the reference.  Converts the 
#information into json format.
#------------------------------------------------------------------------
sub record_reference_info{
	
	my ($reference_id, $species, $genus, $serotype, $sequencer_type, $source, $plasmid_presence, $size, $denovo_assembly, $n50, 
	    $number_of_contigs, $assembly_avg_coverage, $min_contig, $max_contig);
    	        
    my $reference_data = {};
	$reference_data->{'id'} = $reference_id;
	$reference_data->{'species'} = $species;
	$reference_data->{'genus'}= $genus;
	$reference_data->{'serotype'} = $serotype;
	$reference_data->{'sequencer'} = $sequencer_type;
	$reference_data->{'source'} = $source;
	$reference_data->{'assembly'} = $denovo_assembly;
	$reference_data->{'plasmids'} = $plasmid_presence;
	$reference_data->{'size'} = $size;
    $reference_data->{'n50'} = $n50;
    $reference_data->{'number_of_contigs'} = $number_of_contigs;
    $reference_data->{'assembly_avg_coverage'} = $assembly_avg_coverage;
    $reference_data->{'min_contig'} = $min_contig;
    $reference_data->{'max_contig'} = $max_contig;
    
	#if plasmids are present, add their positions to the masked positions file:
    return $reference_data;
    	
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
	my $min_file_size=0;
	my $max_file_size=10;
	my $size_data = {};
    
	foreach(@files){
	    my $output = `ls -l --block-size=M $_`;
	    my $size = (split " ", $output)[4];
	    push @fileSizes, $size;
	    $size_data->{$_} = $size; 	
	}
    #get the rest of hte indivdual file sizes
	foreach(@fileSizes){
		if((split 'M', $_)[0] < $min_file_size){
			#throw warning to user and suggest top up data required
		}
		elsif((split 'M', $_)[0] > $max_file_size){
			#downsample the file 
		}
	}
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

