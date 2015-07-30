package Reporter;

use strict;
use warnings;
use Hash::Merge qw( merge );
use String::Util 'trim';
use JSON;
use Test::JSON;
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
sub bam_quality_data{
    
    my($self, $json_daisy_chain, @bams) = @_;
       
    #the variables to set that will be recorded
    my $reference_size;
    my @warnings;
    my $status = 'PASSED';
    my @common_problems;
    my @common_solutions;
    my %output_data;
    my $depth;
        
    @common_problems = ['Insufficient amount of sequencing data', 'Sample contamination', 'Imporperly identified/labelled samples', 'Poor reference strain used'];                                             
                
    @common_solutions = ['Top up the sequencing data for the poorly mapped isolates', 'Remove isolates that are mapping poorly from the analysis', 'Ensure the samples are identified and labelled correctly', 'Select a different reference strain to map to'];
                             
    #record thresholds used for depth of coverage (default = 15x)
    #parse through the % mapped for each bam file:
    my $bam_string;
	my $x = 1;
    foreach(@bams){
	    $bam_string .= ' --bam bam'.$x.'='.$_;
	    $x++;
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
			chomp $_;
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
    $output_data{'reference'}{'size_from_bam'} = $reference_size;
    $output_data{'bam_stats'}{'min_depth'} = $depth;
    $output_data{'bam_stats'}{'failed_mapping_threshold'} = 80;
    $output_data{'bam_stats'}{'warn_mapping_threshold'} = 90;
    $output_data{'bam_stats'}{'status'} = $status;
    $output_data{'bam_stats'}{'Problem strains'} = \@warnings;
    $output_data{'bam_stats'}{'Common problems'} = \@common_problems;
    $output_data{'bam_stats'}{'Common solutions'} = \@common_solutions;
	
	my $additional = to_json(\%output_data);
    my $output = merge_json($json_daisy_chain, $additional);    
    
    return $output;
}

#------------------------------------------------------------------------
#Purpose:
#    Display the output from the filter-stats.pl script to the user in an
#easily readable format for the reporter log file.  Converts the data into 
#json format as well.
#------------------------------------------------------------------------
sub record_filter_stats{
	
	my($self, $json_daisy_chain, $pseudoalign_filepath) = @_;
	my $result = `perl ../filter-stats.pl -i $pseudoalign_filepath`;
	my %filter_stats;
	
    for(split /^/, $result){
    	switch($_){
    		case {$_ =~ 'generate phylogeny'} {$filter_stats{'filter_stats'}{'total_sites'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'sites identified'} {$filter_stats{'filter_stats'}{'sites_unfiltered'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'sites filtered'} {$filter_stats{'filter_stats'}{'sites_filtered'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'Coverage filtered'} {$filter_stats{'filter_stats'}{'filtered_coverage'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'mpileup filtered'} {$filter_stats{'filter_stats'}{'filtered_mpileup'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'Density filtered'} {$filter_stats{'filter_stats'}{'filtered_density'}=trim((split ":", $_ )[1])}
    		case {$_ =~ 'Invalid filtered'} {$filter_stats{'filter_stats'}{'filtered_invalid'}=trim((split ":", $_ )[1])}
        } 
    }
    
    my $additional = to_json(\%filter_stats);
    my $output = merge_json($json_daisy_chain, $additional);    
    
    return $output;
}

#------------------------------------------------------------------------
#Purpose:
#   Prints user-provided information about the reference.  Converts the 
#information into json format.
#------------------------------------------------------------------------
sub record_reference_info{
	
	my ($self, $json_daisy_chain, $reference_file, $sequencer, $source, $plasmid_presence) = @_;
    
    my $ref_stats = `perl ../ref_stats.pl -i 1000 $reference_file`;
            
    my %reference_data;
    $reference_data{'reference'}{'sequencer'} = $sequencer;
	$reference_data{'reference'}{'source'} = $source; 
	$reference_data{'reference'}{'plasmids'} = $plasmid_presence;
	
	#TODO:
    #$reference_data{'reference'}{'min_contig'} = $min_contig;
    #$reference_data{'reference'}{'max_contig'} = $max_contig;
    
    my $additional = merge_json(to_json(\%reference_data), $ref_stats);
    my $output = merge_json($json_daisy_chain, $additional);    
    
    return $output;
}

#------------------------------------------------------------------------
#Purpose:
#    Record the size of all read/alignment files found in the pipeline/
#------------------------------------------------------------------------
sub record_file_sizes{
		
    my($self, $type, $json_daisy_chain, @files) = @_;
     	
	my @fileSizes;
	my $min_file_size=0;
	my $max_file_size=1000;
	my %size_data;
	
    foreach(@files){
	    my $output = `ls -l --block-size=M $_`;
	    my $size = (split " ", $output)[4];
        push @fileSizes, $size;
        $size_data{'file_sizes'}{$type}{$_} = $size; 	
	}
	
    #get each of the indivdual file sizes
	foreach(@fileSizes){
		if((split 'M', $_)[0] < $min_file_size){
			#throw warning to user and suggest top up data required
		}
		elsif((split 'M', $_)[0] > $max_file_size){
			#downsample the file 
		}
	}
	
	my $additional = to_json(\%size_data);
    my $output = merge_json($json_daisy_chain, $additional);
	
	return $output;
}

#------------------------------------------------------------------------
#Purpose:
#    Record all of the parameters used for the various steps in the 
#pipeline.  All information required to replicate the pipeline run.
#------------------------------------------------------------------------
sub record_run_parameters{	
	
	my($self, $json_daisy_chain, $drmaa_general, $drmaa_trimClean, 
	$drmaa_vcf2core, $drmaa_vcf2pseudoalign, $freebayes_params, 
	$max_coverage, $min_coverage, $mode, $processors, $smalt_index, 
	$smalt_map,	$trim_clean, $vcf2core_cpus, $vcf2pseudo_cpus, $id) = @_;
	
	my %parameters_data;
	$parameters_data{'parameters'}{'id'} = $id;
	$parameters_data{'parameters'}{'drmaa'}{'general'} = $drmaa_general;
	$parameters_data{'parameters'}{'drmaa'}{'trimClean'} = $drmaa_trimClean;
	$parameters_data{'parameters'}{'drmaa'}{'vcf2core'} = $drmaa_vcf2core;
	$parameters_data{'parameters'}{'drmaa'}{'vcf2pseudoalign'} = $drmaa_vcf2pseudoalign;
	$parameters_data{'parameters'}{'freebayes'} = $freebayes_params;
	$parameters_data{'parameters'}{'max_coverage'} = $max_coverage;
	$parameters_data{'parameters'}{'min_coverage'} = $min_coverage;
	$parameters_data{'parameters'}{'mode'} = $mode;
	$parameters_data{'parameters'}{'processors'} = $processors;
	$parameters_data{'parameters'}{'smalt_index'} = $smalt_index;
	$parameters_data{'parameters'}{'smalt_map'} = $smalt_map;
	$parameters_data{'parameters'}{'trim_clean_params'} = $trim_clean;
	$parameters_data{'parameters'}{'vcf2core_cpus'} = $vcf2core_cpus;
	$parameters_data{'parameters'}{'vcf2pseudoalign_cpus'} = $vcf2pseudo_cpus;
	
	my $additional = to_json(\%parameters_data);
    my $output = merge_json($json_daisy_chain, $additional);    
    
    return $output;
}

#------------------------------------------------------------------------
#Purpose:
#        
#------------------------------------------------------------------------
sub vcf2core_stats{
	my($self, $json_daisy_chain, $file_location) = @_;
	
	my %vcf2core_stats;
	
	die "Unable to find the vcf2core.out file at $file_location." if not -e $file_location;
	open(FILE, '<', $file_location) or die "Unable to open the vcf2core.out file.";
	
    while(<FILE>){
    	if(!($_ =~ '#')){
    	   chomp($_);
    	   my @line = split(',', $_);
    	   if( 0+@line == 5 ){
    	       $vcf2core_stats{'vcf2core_stats'}{$line[0]}{'total_length'} = $line[1];
    	       $vcf2core_stats{'vcf2core_stats'}{$line[0]}{'invalid_pos'} = $line[2];
    	       $vcf2core_stats{'vcf2core_stats'}{$line[0]}{'total_core'} = $line[3];
               $vcf2core_stats{'vcf2core_stats'}{$line[0]}{'percent_in_core'} = $line[4];  
    	   }  
    	}
    }
    
	my $additional = to_json(\%vcf2core_stats);
    my $output = merge_json($json_daisy_chain, $additional);    
    
    return $output;			
}

#------------------------------------------------------------------------
#Purpose:
#    Reports the gene within which a particular variant occurs and the 
#type of change that the variant will cause to the nucleotide sequence
#of the resultant protein, if applicable.
#------------------------------------------------------------------------
sub record_snp_eff{
	#IMPLEMENTED TO v2.0
}

#-----------------------------------------------------------------------
#Purpose:
#    Provides a generic plugin interface that allows developers to add 
#additional QC scripts into various stages of the pipeline without
#altering the main Reporter.pm module.
#-----------------------------------------------------------------------
sub run_plugin{
	my($self, $json_daisy_chain, $script_name, $script_variables);
	
	my $command = $script_name;
	
	for(keys in $script_variables){
		#add each command line variable to the command string
		$command = $command.' $_ $script_variables{$_} ';
	}
	
	my %output = `$command`;
	#ensure that the output is in JSON format
	die "Improperly formatted JSON output for $script_name output" if !(is_valid_json(\%output));
	
	my $additional = to_json(\%output);
    my $output = merge_json($json_daisy_chain, $additional);    
    
    return $output;	
}

#------------------------------------------------------------------------
#Purpose:
#    Function to merge two json strings into a single properly formatted
#json string.
#Inout:
#    $original: The first json formatted string.
#    $additional: The second json formatted string.
#Output:
#    A JSON string representing the union of the two original json
#strings.
#------------------------------------------------------------------------
sub merge_json{
	my($original, $additional) = @_;
    
    my $first = from_json($original);
    my $second = from_json($additional);
    
    my %output = %{merge($first, $second)};
    
    return to_json(\%output, { pretty => 1 }); 
}

1;