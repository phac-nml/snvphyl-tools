#!/usr/bin/env perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use autodie;
use FindBin qw($Bin);
use File::Basename;
use Parallel::ForkManager;
use Readonly;

Readonly my $MIN_DEPTH => 10;    # THIRD COLUMN
Readonly my $MIN_MAP => 80;

__PACKAGE__->run unless caller;


sub run {
    my ( $size, $man, $help, $min_depth, $log_dir, $min_map, %bam_files, $cores );
	
    GetOptions(
        "c|cores=i"   => \$cores,
        "bam=s"		  => \%bam_files,
        "l|log_dir=s" => \$log_dir,
        "s|size=s"    => \$size,
        "min-map=f"	  => \$min_map,
        "min-depth=i" => \$min_depth,
        "h|help"      => \$help,
        "m|man"       => \$man
    );
    pod2usage(1) if $help;
    pod2usage(-verbose => 2) if $man;
	
    unless ( (scalar keys %bam_files != 0 ) ) {
        print "Unable to find any input bam files.\n\n";
        pod2usage(1);
    }
    unless ( defined $size ) {
        print "Calculating genome size from input BAM files.\n\n";
    }
		
	#set default values if undefined on command line
    if(!defined $min_depth){
        $min_depth = $MIN_DEPTH;
        print "  Using a minumum depth of $MIN_DEPTH\n";
        print "  Use the --min-depth flag to set custom depth\n\n";
    }
    else{
    	print "Using a min_depth of ".$min_depth."\n";
    }
    #set default number of cores
	$cores = 1 if (not defined $cores);
	#set default log directory to local directory    
	$log_dir="./" if (not defined $log_dir);
	#set default minimum percent mapping
	$min_map=$MIN_MAP if (not defined $min_map);
	
	#create the log file to print warnings to
	open(my $log, '>', $log_dir.'mapping_percentage.log');
	
	#retrieve all of the bam file locations from the hash
    my @files = values %bam_files;
 
 	#ensure that bam files are properly input on the command line and that each file path exists
 	if (@files <= 0){ die "Error: No bam files input."};
 	foreach(@files){
 		if (!-e $_) {die "Error: Invalid bam file referenced."};
 	}
 	
	#command to get the size of the genome from the bam file:
	$size = `samtools view -H $files[0] | grep -P '^\@SQ' | cut -f 3 -d ':' | awk '{sum+=\$1} END {print sum}'`;
	
	die "Error: Size of reference genome could not be determined." if (not defined $size);	
	#ensure the size of the genome is properly determined
	
	#parse results and determine what should be written for user to view
	my @results;
	@results = verify_percent_coverage( \@files, $size, $min_depth, $cores );
	print $log "==========Reference Mapping Quality===========\n";
	print $log "NUMBER OF BP's IN REFERENCE GENOME: ".$size."\n";
    foreach my $result(@results){
    	my @split = split(',', $result);
    	my @double = split('%', $split[1]);
    	print $log "Mapping to reference for isolate ".$split[0]." is ".$split[1]."\n" if $double[0] < $min_map; 
    }
	   	
}

#----------------------------------------------------------#
# Calculate the % Coverage for all .bam files in directory #
#----------------------------------------------------------#
sub verify_percent_coverage {
    my ( $files, $size, $min_depth, $cores ) = @_;

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
		
        # Now that we have the results...
        #--------------------------------#
        my $gap_length = get_gap_length($result, $min_depth);

        my $line;
        if ( $gap_length == -1) {
            $line = "$name,0%";
        }
        else {
            $line = sprintf "$name,%3.2f%%", (( $size - $gap_length ) / $size * 100);            
        }


        $pm->finish(0,\$line);
    }
    $pm->wait_all_children;
	
    #return the results
    return @results;
}

#----------------------------------#
# Find the length with < MIN_DEPTH #
#----------------------------------#
sub get_gap_length {
    my ($result, $min_depth) = @_;

    my @lines = split /\n/, $result;

    my $previous_pos = 0;
    my $gap_total    = -1;
    foreach (@lines) {
        my ( undef, $pos, $count ) = split /\s+/;

        # OOPS, You skipped an area
        #--------------------------#
        if ( $previous_pos != $pos - 1 ) {
            my $start = $previous_pos + 1;
            my $end   = $pos - 1;
            $gap_total += $end - $start;

        }
        elsif ( $count <= $min_depth ) {
            $gap_total++;
        }
        $previous_pos = $pos;
    }
    return $gap_total;
}

=pod

=head1 NAME

verify_mapping_quality.pl - Script to check the mapping quality of all BAM files generated in the core SNP pipeline.

=head1 VERSION

This documentation refers to verify_mapping_quality.pl version 0.0.1.

=head1 SYNOPSIS

verify_mapping_quality.pl -l /log-direcotry -i /inputDataDirectory --min-depth minimum-depth --min-map minimum-percent-mapping -s genome-size(bp) -h help

=head1 OPTIONS

=over

=item B<-l>, B<--log-dir> [required]

The directory location where the log file exists.

=item B<-i>, B<--input-dir> [optional]

The input directory that contains all of the BAM files for the run.

=item B<--min-depth> [optional]

The minimum depth of coverage required in each BAM file, pipeline will die unless this minimum depth is met.

=item B<--min-map> [optional]

The minimum percent coverage required in each BAM file, pipeline will log all BAM files that do not meet the coverage limit.

=item B<-s>, B<--genome-size> [optional]

The minimum percent coverage required in each BAM file, pipeline will die unless this minimum percent mapping is met.

=item B<-h>, B<--help>

To display help message

=back

=head1 DESCRIPTION

verify_mapping_quality will check BAM files to check the percent of the isolate genomes that have been mapped to the reference sequence. 

=cut
