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

Readonly my $MIN_DEPTH => 10;    
Readonly my $MIN_MAP => 80;

__PACKAGE__->run unless caller;


sub run {
    my ( $size, $man, $help, $min_depth, $min_map, %bam_files, $cores ,$output);
	
    GetOptions(
        "c|cores=i"   => \$cores,
        "bam=s"		  => \%bam_files,
        "s|size=s"    => \$size,
        "min-map=f"	  => \$min_map,
        "min-depth=i" => \$min_depth,
        'output=s'       => \$output,
        "h|help"      => \$help,
        "m|man"       => \$man
    );
    pod2usage(0) if $help;
    pod2usage(-verbose => 2) if $man;
	
    unless ( (scalar keys %bam_files != 0 ) ) {
        print STDERR "Unable to find any input bam files.\n\n";
        pod2usage(0);
    }

	#set default values if undefined on command line
    if(!defined $min_depth){
        $min_depth = $MIN_DEPTH;
    }
			
    #set default number of cores
    $cores = 1 if (not defined $cores);
    #set default minimum percent mapping
    $min_map=$MIN_MAP if (not defined $min_map);

    #check to see if we are given a $out , if not we will write to STDOUT
    my $out_fh;
    if ( defined $output ) {
        open( $out_fh, '>', $output );
    }
    else {
        $out_fh = \*STDOUT;
    }
    
    
    #retrieve all of the bam file locations from the hash
    my @files = values %bam_files;
 
 	#ensure that bam files are properly input on the command line and that each file path exists
 	if (@files <= 0){ die "Error: No bam files input."};
 	foreach(@files){
 		if (!-e $_) {die "Error: Invalid bam file referenced."};
 	}
 	
	#command to get the size of the genome from the bam file:
	my @sizeArray = `samtools view -H $files[0] | grep -P '^\@SQ' | cut -f 3 -d ':'`;
	
	#if the read maps to several reference contigs, then add the lengths of each reference
	#contig to calculate the total length for the reference	
	$size = 0;
	foreach(@sizeArray){ 
	   $size += $_;
	}	
	die "Error: Size of reference genome could not be determined." if (not defined $size || $size eq 0);	
		
	#parse results and determine what should be written for user to view
	my @results;

	@results = verify_percent_coverage( \%bam_files, $size, $min_depth, $cores );
	print $out_fh "==========Reference Mapping Quality===========\n";
	print $out_fh "NUMBER OF BP's IN REFERENCE GENOME: ".$size."\n";
	print $out_fh "MINIMUM DEPTH: ".$min_depth."\n";
	print $out_fh "MINIMUM MAPPING: ".$min_map."\n";
	
    foreach my $result(@results){
    	my @split = split(',', $result);
    	my @double = split('%', $split[1]);
    	print $out_fh $split[0]." : ".$split[1]."\n" if $double[0] < $min_map; 
    }

    return;
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
    
    foreach my $file ( keys %$files ) {
        $pm->start && next;

        my $name = $file;

        # Run samtools depth and get results
        #------------------------------------#
        my $input_file = $files->{$file};
        my $result = `samtools depth $input_file`;
        
        #check for errors that occur while running samtools depth
		die "Error: samtools depth exited with error while working with $file.\n" if (not defined $result);
		
        # Now that we have the results...
        #--------------------------------#
        my $total_passed = total_passed_positions($result, $min_depth);
		
        my $line = sprintf "$name,%3.2f%%", (($total_passed/ $size) * 100);  

        $pm->finish(0,\$line);
    }
    $pm->wait_all_children;
	
    #return the results
    return @results;
}

#
#Returns the total number of positions that pass the min depth
#
sub total_passed_positions{
	my ($result, $min_depth) = @_;
    my $total_passed = 0;
    
    my @lines = split /\n/, $result;
    foreach (@lines){
    	my ( undef, $pos, $count ) = split /\s+/;
    	if($count >= $min_depth){
    		$total_passed++;
    	}
    }
    return $total_passed;    
}


=pod

=head1 NAME

verify_mapping_quality.pl - Script to check the mapping quality of all BAM files generated in the core SNP pipeline.

=head1 VERSION

This documentation refers to verify_mapping_quality.pl version 0.0.1.

=head1 SYNOPSIS

verify_mapping_quality.pl --bam bamX=/inputDirrectory/bamfile.bam --min-depth minimum-depth --min-map minimum-percent-mapping -h help

=head1 OPTIONS

=over

=item B<--bam> [REQUIRED]

The location for a specific BAM file in the dataset. Multiple BAM files can be input.  Example with 3 BAM files: --bam bam1=/path/bam1.bam --bam bam2=/path/bam2.bam --bam bam3=/path/bam3.bam

=item B<--min-depth> [optional]

The minimum depth of coverage required at each genome position to be considered mapped.  Default value is 15x.

=item B<--min-map> [optional]

The minimum percent mapped to reference for each strain, pipeline will log all strains that do not meet this minimum percentage. Default value is 80%.

=item B<-h>, B<--help>

To displays help screen.

=back

=head1 DESCRIPTION

verify_mapping_quality will check BAM files and list any samples that map to <= --min-map of the reference sequence. 

=cut
