use FindBin;
use strict;
use warnings;
use Hash::Merge qw( merge );
use String::Util 'trim';
use JSON;
use Test::JSON;
use Switch;
use Getopt::Long;

my $script_dir = $FindBin::Bin;
our $VERSION = 0.1;

__PACKAGE__->run unless caller;

sub run{
    
   my ( $step, $output_json, $json_input, %bam_files, $pseudoalign_filepath, $reference_filepath, $ref_sequencer, $ref_source, 
        $plasmids, $genus, $species, $serotype, $file_type, %file_sizes, $drmaa_general, $drmaa_trimClean, 
        $drmaa_vcf2core, $drmaa_vcf2pseudoalign, $freebayes_params, $max_coverage, $min_coverage, $mode, 
        $processors, $smalt_index, $smalt_map, $trim_clean, $vcf2core_cpus, $vcf2pseudo_cpus, $id, $masked_positions, 
        %read_files, $vcf2core_stats);
 
   GetOptions(
      "step=s" => \$step,
      "json=s" => \$json_input,
      "output=s" => \$output_json,
      "bam=s" => \%bam_files,
      "pseudo=s" => \$pseudoalign_filepath,
      "ref-file=s" => \$reference_filepath,
      "ref-sequencer=s" => \$ref_sequencer,
      "ref-source=s" => \$ref_source,
      "plasmids=s" => \$plasmids,
      "genus=s" => \$genus,
      "species=s" => \$species,
      "serotype=s" => \$serotype,
      "file-type=s" => \$file_type,
      "file-sizes=s" => \%file_sizes,
      "freebayes-params=s" => \$freebayes_params,
      "max-coverage=i" => \$max_coverage,
      "min-coverage=i" => \$min_coverage,
      "mode=s" => \$mode,
      "processors=i" => \$processors,
      "smalt-index=s" => \$smalt_index,
      "smalt-map=s" => \$smalt_map,
      "trim-clean=s" => \$trim_clean,
      "vcf2core-cpus=i" => \$vcf2core_cpus,
      "run-id=s" => \$id,
      "masked-positions=s" => \$masked_positions,
      "read-file=s" => \%read_files,
      "vcf2core-stats=s" => \$vcf2core_stats
   );
   
   my $json_daisy_chain="";
   my $output;
   
   if(defined $json_input){
      open(JSON_FILE, $json_input) or die "Unable to open input file handle.";  
   	  while(<JSON_FILE>){
   	     $json_daisy_chain .= $_;
   	     print $_;	
   	  }
   }
   else{
      die "Input json file, $json_input, not defined.";
   }
   
   switch($step){
         case "bam_quality_data" { $output = bam_quality_data($json_daisy_chain, %bam_files) }
         case "record_filter_stats"{ $output = record_filter_stats($json_daisy_chain, $pseudoalign_filepath)}
         case "record_reference_info"{ $output = record_reference_info($json_daisy_chain, $reference_filepath, $ref_sequencer, $ref_source, $plasmids, $genus, $species, $serotype)};
         case "record_file_sizes"{ $output = record_file_sizes($json_daisy_chain, $file_type, %file_sizes)};
         case "record_run_parameters"{ $output = record_run_parameters($json_daisy_chain, $freebayes_params, $max_coverage, $min_coverage, $mode, $processors, $smalt_index, $smalt_map, $trim_clean, $vcf2core_cpus, $id, 
                                                    $masked_positions, %read_files)};
         case "vcf2core_stats"{$output = vcf2core_stats($json_daisy_chain, $vcf2core_stats)};
         #TODO:case {'run_plugin'}{return run_plugin($json_daisy_chain)};
   }
   #print the json output to a file for the results
   open(my $out, '>', $output_json);
   print $out $output;
   close(JSON_FILE);
   close($out);
}

#------------------------------------------------------------------------
#Purpose:
#    Analyzes the quality with which each .bam file mapped to the 
#reference sequence for the run.  Marks the results as PASS, WARN, or 
#FAIL depending on the % mapped for each .bam input.  
#Input:
#   $json_daisy_chain: JSON string containing reporter data for the run.
#   @bams: A list of file paths to .bam files.
#Output:
#   JSON daisy chain with bam quality data included. 
#------------------------------------------------------------------------
sub bam_quality_data{
    
    my($json_daisy_chain, %bams) = @_;

    #the variables to set that will be recorded
    my $reference_size;
    my @warnings;
    my $status = "PASSED";
    my @common_problems;
    my @common_solutions;
    my %output_data;
    my $depth;
        
    @common_problems = ['Insufficient amount of sequencing data', 'Sample contamination', 'Imporperly identified/labelled samples', 'Poor reference strain used'];                                             
                
    @common_solutions = ['Top up the sequencing data for the poorly mapped isolates', 'Remove isolates that are mapping poorly from the analysis', 'Ensure the samples are identified and labelled correctly', 'Select a different reference strain to map to'];
                             
    #record thresholds used for depth of coverage (default = 15x)
    #parse through the % mapped for each bam file:
    my $bam_string;
    
    for(keys %bams){
        $bam_string .= ' --bam '.$_.'='.$bams{$_};
    }
    
    my $mapping_results = `perl $script_dir/../verify_mapping_quality.pl $bam_string --min-map 90 --min-depth 15`;
               
    #parse the results from the verify_mapping_quality.pl script
    for(split /^/, $mapping_results){    
        if($_ =~ 'MINIMUM DEPTH'){
            $depth = trim((split ':', $_)[1]);
            chomp $depth;
        }
        if($_ =~ 'REFERENCE GENOME'){
            $reference_size = trim((split ':', $_)[1]);
            chomp $reference_size;
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
            $status = "FAILED";
        }
        elsif(($percent_mapped <= 90) && ($status eq 'PASSED')){
            $status = "WARNING";
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
#    Analyze the output from the filter-stats.pl script to the user in an
#easily readable format for the reporter log file.
#Input:
#   $json_daisy_chain: JSON string containing reporter data for the run.
#   $pseudoalign_filepath: File path to the pseudoalign-positions.tsv
#                          for the run.
#Output:
#   JSON daisy chain with the filter stats data included.
#------------------------------------------------------------------------
sub record_filter_stats{
    
    my($json_daisy_chain, $pseudoalign_filepath) = @_;
    my $file = $script_dir.$pseudoalign_filepath;
    my $result = `perl $script_dir/../filter-stats.pl -i $file`;
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
#   Records information about the reference used in the run.
#Input:
#   $json_daisy_chain: JSON string containing reporter data for the run.
#   $reference_file: The file location for the reference .fasta file.
#   $sequencer: The type of sequencer used to generate the reference reads.
#   $source: The lab/site/team that the reference sequence was retrieved
#            from.
#   $plasmid_presence: YES/NO flag indicating whether a plasmids are 
#                      present in the reference sequence.
#   $genus: The genus of the reference organism.
#   $species: The species of the reference organism.
#   $serotype: The serotype of the reference organism.
#Output:
#   JSON daisy chain with the reference information included.
#------------------------------------------------------------------------
sub record_reference_info{
    
    my ($json_daisy_chain, $reference_file, $sequencer, $source, $plasmid_presence, $genus, $species, $serotype) = @_;
    
    my $file = $script_dir.'/'.$reference_file;
    my $ref_stats = `perl $script_dir/../ref_stats.pl -i 1000 $file`;
            
    my %reference_data;
    $reference_data{'reference'}{'sequencer'} = $sequencer;
    $reference_data{'reference'}{'source'} = $source; 
    $reference_data{'reference'}{'plasmids'} = $plasmid_presence;
    $reference_data{'reference'}{'genus'} = $genus;
    $reference_data{'reference'}{'species'} = $species;
    $reference_data{'reference'}{'serotype'} = $serotype;
    $reference_data{'input_files'}{'input_reference'} = $reference_file;
    
    #TODO:
    #$reference_data{'reference'}{'min_contig'} = $min_contig;
    #$reference_data{'reference'}{'max_contig'} = $max_contig;
    
    my $additional = merge_json(to_json(\%reference_data), $ref_stats);
    my $output = merge_json($json_daisy_chain, $additional);    
    
    return $output;
}

#------------------------------------------------------------------------
#Purpose:
#    Records information about the size of files used in the pipeline.
#Input:
#   $json_daisy_chain: JSON string containing reporter data for the run.
#   $type: The file format.
#   @files: The list of file paths for files to be analyzed relative to 
#           current working directory.
#Output:
#   JSON daisy chain with the file size information included.   
#------------------------------------------------------------------------
sub record_file_sizes{
        
    my($json_daisy_chain, $type, %files) = @_;
         
    my @fileSizes;
    my $min_file_size=0;
    my $max_file_size=1000;
    my %size_data;
 
    for(keys %files){
        my $file = $script_dir.'/'.$files{$_};
        my $output = `ls -l --block-size=M $file`;
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
#    Records information on all of the parameters used for the various 
#stages of the pipeline. 
#Input:
#   $json_daisy_chain: JSON string containing reporter data for the run.
#   $drmaa_general: drmaa general parameter string.
#   $drmaa_trimClean: drmaa_trimclean parameter string.
#   $drmaa_vcf2core: drmaa_vcf2core parameter string.
#   $drmaa_vcf2pseudoalign: drmaa_vcf2pseudoalign parameter string.
#   $freebayes_params: The freebayes parameter string.
#   $max_coverage: The max coverage value.
#   $min_coverage: The min coverage value.
#   $mode: The mode used for the run.
#   $processors: The number of processors used for parallelization.
#   $smalt_index: The parameter string for smalt_index.
#   $smalt_map: The parameter string for smalt_map.
#   $trim_clean: The parameter string for trim_clean.
#   $vcf2core_cpus: The number of cpus used for the vcf2core stage.
#   $vcf2pseudo_cpus: The number of cpus used for the vcf2pseudo stage.
#   $id: The run id.
#   $masked_positions: The file path to the masked positions file.
#   @read_files: A list of file paths to the fastq files used in the run.
#Output:
#   JSON daisy chain with the run parameter information included.   
#------------------------------------------------------------------------
sub record_run_parameters{    
    
    my($json_daisy_chain, $freebayes_params, 
    $max_coverage, $min_coverage, $mode, $processors, $smalt_index, 
    $smalt_map,    $trim_clean, $vcf2core_cpus, $id,
    $masked_positions, %read_files) = @_;
    
    my %parameters_data;
    $parameters_data{'parameters'}{'id'} = $id;
    $parameters_data{'parameters'}{'freebayes'} = $freebayes_params;
    $parameters_data{'parameters'}{'max_coverage'} = $max_coverage;
    $parameters_data{'parameters'}{'min_coverage'} = $min_coverage;
    $parameters_data{'parameters'}{'mode'} = $mode;
    $parameters_data{'parameters'}{'processors'} = $processors;
    $parameters_data{'parameters'}{'smalt_index'} = $smalt_index;
    $parameters_data{'parameters'}{'smalt_map'} = $smalt_map;
    $parameters_data{'parameters'}{'trim_clean_params'} = $trim_clean;
    $parameters_data{'parameters'}{'vcf2core_cpus'} = $vcf2core_cpus;
    for(keys %read_files){
        $parameters_data{'input_files'}{'input_bams'}{$_} = $read_files{$_}; 
    } 
    $parameters_data{'input_files'}{'masked_positions'} = $masked_positions;
    
    my $additional = to_json(\%parameters_data);
    my $output = merge_json($json_daisy_chain, $additional);    
    
    return $output;
}

#------------------------------------------------------------------------
#Purpose:
#   Script to analyze and report the vcf2core stats from the pipeline
#run.
#Input:
#   $json_daisy_chain: JSON string containing reporter data for the run.
#   $file_location: The file location for the vcf2core stats output file.
#Output:
#   JSON daisy chain with the vcf2core stats included.   
#------------------------------------------------------------------------
sub vcf2core_stats{
    my($json_daisy_chain, $file_location) = @_;
    
    my %vcf2core_stats;
    
    die "Unable to find the vcf2core.out file at $file_location." if not -e $script_dir.'/'.$file_location;
    open(FILE, '<', $script_dir.'/'.$file_location) or die "Unable to open the vcf2core.out file.";
    
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
    close(FILE);
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
#Input:
#   $json_daisy_chain: JSON string containing reporter data for the run.
#   $script_name: The filepath and name of the script to run.
#   $script_variables: A hash reference that contains a list of script
#                      variables.  Must be in the format of KEY:VALUE,
#                      where KEY = '-x OR --XXXXX' and VALUE= 'XXXXX'
#Output:
#   JSON daisy chain with the plugin information included.   
#-----------------------------------------------------------------------
sub run_plugin{
    my($json_daisy_chain, $plugin_name, $script_variables);
    
    my $command = $plugin_name;
    
    for(keys in $script_variables){
        #add each command line variable to the command string
        $command = $command.' $_ $script_variables{$_} ';
    }
    
    my %output = `$command`;
    #ensure that the output is in JSON format
    die "Improperly formatted JSON output for $plugin_name output" if !(is_valid_json(\%output));
    
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
	
	my %output;
	my $first;
	my $second;
	
	if(!($original eq "")){
       $first = from_json($original);
       $second = from_json($additional);
       %output = %{merge($first, $second)};
       return to_json(\%output);
	}
    else{
    	return $additional;
    }
}

1;