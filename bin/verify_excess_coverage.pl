#!/usr/bin/env perl

use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;
use autodie;
use FindBin qw($Bin);
use File::Basename;
use Parallel::ForkManager;

__PACKAGE__->run unless caller;

#========================================================================
#Purpose:
#    Script to identify regions of the genome that are deemed to have 
#an excessive amount of coverage (peak) based on their std deviation from 
#the mean depth of coverage across all isolates.
#
#Input:
#    $max_stddev -> The maximum +/- number of standard devitaions. 
#                   allowed for a 'normal' depth of coverage.
#    %bam_files  -> Hash of all input bam files.
#    $cores      -> Number of cores to use for paralellization.
#Output:
#    A tsv file describing the areas of excessive or lacking coverage
#for each strain.  The format of each line is as follows:
#    #$strain
#    $contig.\t.$start.\t.$end.\n
#Where:
#    $strain  -> The strain ID.
#    $contig  -> The contig ID.
#    $start   -> The start position for the region.
#    $end     -> The end position for the region.
#========================================================================
sub run{
        
    my ( $man, $help, $max_stddev, $out_dir, %bam_files, $cores );
    
    GetOptions(
        "c|cores=i"   => \$cores,
        "bam=s"           => \%bam_files,
        "max-dev=i" => \$max_stddev,
        "h|help"      => \$help,
        "m|man"       => \$man
    );
    pod2usage(0) if $help;
    pod2usage(-verbose => 2) if $man;
        
    unless ( (scalar keys %bam_files != 0 ) ) {
        print "Unable to find any input bam files.\n\n";
        pod2usage(0);
    }

	#set default values if undefined on command line
    if(!defined $max_stddev){
 		$max_stddev = 5;
    }
    else{
    	print "Using a max std deviation of ".$max_stddev."\n";
    }
    
    #set default number of cores
	$cores = 1 if (not defined $cores);
	
	#retrieve all of the bam file locations from the hash
    my @files = values %bam_files;
 	
 	#ensure that bam files are properly input on the command line and that each file path exists
 	if (@files <= 0){ die "Error: No bam files input."};
 	foreach(@files){
 		if (!-e $_) {die "Error: Invalid bam file referenced at $_."};
 	}
 	
	#command to get the size of the genome from the bam file:
	my @sizeArray = `samtools view -H $files[0] | grep -P '^\@SQ' | cut -f 3 -d ':'`;
	
	#if the read maps to several reference contigs, then add the lengths of each reference
	#contig to calculate the total length for the reference	
	my $size = 0;
	foreach(@sizeArray){ 
	   $size += $_;
	}	
	die "Error: Size of reference genome could not be determined." if (not defined $size || $size eq 0);	
		
	#parse results and determine what should be written for user to view
	my @results;
	
	#retrieve the depths of coverage for all positions in each of the strains
	@results = depth_of_coverage( \%bam_files, $size, $cores, $max_stddev );    
    
    print @results;
    
    exit;
}

sub calc_mean_coverage{
    my(@data) = @_;
    if (not @data) {
        die("Empty array.");
    }
    my $count = 0;
    my $total = 0;
    for(split /^/, $data[0]) {
        my @tab_parsed = split('\t', $_);
        if($tab_parsed[2] && ($tab_parsed[2] >= 0)){
            $total += $tab_parsed[2];
            $count++;
        }
    }
    my $average = $total / $count;
    return $average;
}

sub calc_stddev{
	my($average, @data) = @_;

	my $count = 0;
    my $sqtotal = 0;
    for(split /^/, $data[0]) {
        my @tab_parsed = split('\t', $_);
        if($tab_parsed[2] && ($tab_parsed[2]>=0)){
            $sqtotal += ($average - $tab_parsed[2]) ** 2;
            $count ++;
        }
    }
    my $std = ($sqtotal / ($count - 1)) ** 0.5;
    return $std;
}

#========================================================================
#Purpose:
#    Algorithm to identify excesive coverage regions within a genome. 
#Input:
#    $max_dev -> Any position with over max deviation coverage depth is 
#                considered an excessive coverage position.
#    @rows     -> List of genome regions with excessive coverage
#Output:
#    $result -> A String describing a list of high density regions in 
#               the genome.
#========================================================================
sub find_excessive_coverage{
	my ($threshold, @rows) = @_;
	
	my $result = "";
	my $current_chromosome = '';
	my $density_flag = 0;
	my $start = 0;
	my $end = 0;
	my $depth_of_coverage=0;
	
	foreach(@rows){
		for(split /^/, $_){
		my @parsed_row = split '\t', $_;
		#if the chromosome/contig changes, do the same as above
		if(!($current_chromosome eq $parsed_row[0])){
			if($density_flag){
			    $result .= "$current_chromosome\t$start\t$end\n";
			}
			$current_chromosome = $parsed_row[0];
	        $density_flag = 0;
	        $start = $parsed_row[1];
	        $end = $parsed_row[1];
		}
		#check if the current position has excessive coverage
		if((exists $parsed_row[2]) && ($parsed_row[2] > $threshold)){
			$density_flag = 1;
		}
		elsif(exists $parsed_row[2]){
			if($density_flag){
			    $result .= "$parsed_row[0]\t$start\t$end\n";
			}
			$start = $parsed_row[1];
			$density_flag = 0;
		}
		$end = $parsed_row[1];
	}
	}
	return $result;
}

#========================================================================
#Purpose:
#    Method to retrieve the depth of coverage for each position in every
#input .bam alignment file and then find regions of excessive coverage 
#for each isolate. 
#Input:
#    $files -> An array of bam file locations
#    $size  -> The size of the reference genome
#    $max_stddev -> The maximum # of std deviations from mean that 
#                   'normal' coverage can be.
#    $cores -> The number of cores to be used for parallelization.
#Output:
#    A list of tsv lines describing each region of excessive coverage.
#Format for each line is: $contig.\t.$start.\t.$end\n. Strains 
#are indicated by a newline in the form: #$strain\n.
#========================================================================
sub depth_of_coverage {
    my ( $files, $size, $cores, $max_stddev ) = @_;

    # Create a Prallel::ForkManager
    #------------------------------#
    my $pm = new Parallel::ForkManager($cores);

    my @results;
    
    # data structure retrieval and handling
    $pm -> run_on_finish ( # called BEFORE the first call to start()
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
            # retrieve data structure from child
            if (defined($data_structure_reference)) {  # children are not forced to send anything
                push @results,${$data_structure_reference};
            } else {  # problems occuring during storage or retrieval will throw a warning
                die "Did not recieve data from one of the children\n";
            }
        }
    );
    
    # Let's do the heavy lifing
    #--------------------------#
    foreach my $file ( keys %$files ) {
        $pm->start && next;

		my $name = $file;

        # Run samtools depth and get results
        #------------------------------------#
        my $input_file = $files->{$file};
        my $result = `samtools depth $input_file`;
        
        #check for errors that occur while running samtools depth
		die "Error: samtools depth exited with error while working with $input_file.\n" if (not defined $result);
		
		#calculate the averages and std deviation for the coverage stats of each strain
        my $average = calc_mean_coverage($result);
        my $std_dev = calc_stddev($average, $result);
        #print to two decimal places:
        $std_dev = sprintf "%.2f", $std_dev;
        my $max_dev = sprintf "%.0f", ($average + ($max_stddev * $std_dev));
        
		my $excessive_coverage = find_excessive_coverage($max_dev, $result);
		#add the strain name as the first line in the set of results
		my $strain = (split('\.', $name))[0];
		$result = "#$strain\n".$excessive_coverage;
		
        $pm->finish(0,\$result);
    }
    $pm->wait_all_children;
	
    #return the results
    return @results;
}

=pod

=head1 NAME

verify_excess_coverage.pl - Script to check for regions of excess coverage within isolate genomes.  

=head1 SYNOPSIS

verify_excess_coverage.pl -c [NUM_CPU] --max-dev [MAX_STD_DEV] --bam bam1=/path/to/bam1 --bam bam2=/path/to/bam2 --bam bamX=/path/to/bamX

=head1 OPTIONS

=over

=item B<--bam> [REQUIRED]

The location for a specific BAM file in the dataset. Multiple BAM files can be input.  Example with 3 BAM files: 
	
--bam bam1=/path/bam1.bam --bam bam2=/path/bam2.bam --bam bam3=/path/bam3.bam

=item B<-c>, B<--cores> [optional]                                                                                                      

The number of CPU cores that should be used for the calculations.

=item B<--max-dev> [optional]

The maximum number of standard deviations from the mean allowable for positions with a 'normal' depth of coverage.

=item B<-h>, B<--help>

To displays help screen.

=back

=head1 DESCRIPTION

verify_excess_coverage - Script to identify regions of the genome that are deemed to have an excessive depth of coverage (peak) based on their std deviation from the mean depth.

=cut