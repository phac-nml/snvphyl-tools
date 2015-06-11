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
#    Script to identify positions in isolate alignments that have low or
#zero coverage.
#
#Input:
#    $max_depth  -> The maximum depth of coverage for low coverage po
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
        
    my ( $man, $help, $max_depth, %bam_files, $cores );
    
    GetOptions(
        "c|cores=i"   => \$cores,
        "bam=s"           => \%bam_files,
        "max-depth=i" => \$max_depth,
        "h|help"      => \$help,
        "m|man"       => \$man
    );
    pod2usage(1) if $help;
    pod2usage(-verbose => 2) if $man;
        
    unless ( (scalar keys %bam_files != 0 ) ) {
        print "Unable to find any input bam files.\n\n";
        pod2usage(1);
    }

	#set default values if undefined on command line
    if(!defined $max_depth){
 		$max_depth = 0;
    }
    else{
    	print "Using a max depth of ".$max_depth."\n";
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
	
	#retrieve the depths of coverage for all positions in each of the strains
	my @results = check_low_coverage( \@files, $cores, $max_depth );    
    
    print @results;
    
    exit;
}


#========================================================================
#Purpose:
#    Function identifies regions in isolate genomes that have low depth
#coverage as defined by the user.
#Input:
#    $files -> An array of bam file locations
#    $max_depth -> The maximum depth for positions considered 'low depth'
#    $cores -> The number of cores to be used for parallelization.
#Output:
#    A list of tsv lines describing each region of excessive coverage.
#Format for each line is: $contig.\t.$start.\t.$end\n. Strains 
#are indicated by a newline in the form: #$strain\n.
#========================================================================
sub check_low_coverage {
    my ( $files, $cores, $max_depth ) = @_;

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
    foreach my $file ( @$files ) {
        $pm->start && next;
        
        # Get the actual file name
        #-------------------------#
        my @suffix = [".bam", ".dat"];
        my $name = fileparse( $file, @suffix );

        # Run samtools depth and get results
        #------------------------------------#
        my $result = `samtools depth $file`;
        
        #check for errors that occur while running samtools depth
		die "Error: samtools depth exited with error while working with $file.\n" if (not defined $result);
        
        my $low_coverage;
        my $low_coverage_flag = 0;
        my $low_coverage_start;
        my $low_coverage_end;
        foreach($result){
        	for(split /^/, $_){
        	my @tab_parsed = split('\t', $_);
        	if($tab_parsed[2] <= $max_depth){
        		if(!$low_coverage_flag){
        			$low_coverage_start = $tab_parsed[1];
        		}
        		$low_coverage_end = $tab_parsed[1];
        		$low_coverage_flag = 1;
        	}
        	else{
        		if($low_coverage_flag){
        			$low_coverage .= "$tab_parsed[0]\t".$low_coverage_start."\t".$low_coverage_end."\n";
        		}
        		$low_coverage_flag = 0;
        	}
        }
        }
		#add the strain name as the first line in the set of results
		my $strain = (split('\.', $name))[0];
		$result = "#$strain\n".$low_coverage;
		
        $pm->finish(0,\$result);
    }
    $pm->wait_all_children;
	
    #return the results
    return @results;
}

=pod

=head1 NAME

verify_zero_depth.pl - Script to check for regions of low coverage depth within isolate genomes.  

=head1 SYNOPSIS

verify_low_depth.pl -c [NUM_CPU] --max-depth [MAX_DEPTH] --bam bam1=/path/to/bam1 --bam bam2=/path/to/bam2 --bam bamX=/path/to/bamX

=head1 OPTIONS

=over

=item B<--bam> [REQUIRED]

The location for a specific BAM file in the dataset. Multiple BAM files can be input.  Example with 3 BAM files: 
	
--bam bam1=/path/bam1.bam --bam bam2=/path/bam2.bam --bam bam3=/path/bam3.bam

=item B<-c>, B<--cores> [optional]                                                                                                      

The number of CPU cores that should be used for the calculations.

=item B<--max-depth> [optional]

The maximum depth for a position that is considered to have a 'low' depth of coverage.

=item B<-h>, B<--help>

To displays help screen.

=back

=head1 DESCRIPTION

verify_low_depth - Script to identify regions of the isolate sequences that are deemed to have a low depth of coverage (valley).

=cut	