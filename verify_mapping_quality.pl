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

1;

sub run {
    my ( $dir, $size, $man, $help, $min_depth, $log_dir, $min_map );

    GetOptions(
        "i|dir=s"     => \$dir,
        "l|log_dir=s" => \$log_dir,
        "s|size=s"    => \$size,
        "min-map=f"	  => \$min_map,
        "min-depth=i" => \$min_depth,
        "h|help"      => \$help,
        "m|man"       => \$man
    );
    pod2usage(1) if $help;
    pod2usage(-verbose => 2) if $man;

    unless ( (defined $dir && -e $dir) ) {
        print "Please specifiy a directory or file name.\n\n";
        pod2usage(1);
    }
    unless ( defined $size ) {
        print "Please specify the size of the genome.\n\n";
        pod2usage(1);
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
	$log_dir='' if (not defined $log_dir);
	$min_map=$MIN_MAP if (not defined $min_map);
	
	#create the log file to print warnings to
	open(my $log, '>'.$log_dir.'mapping_percentage.log');
	
    my @files;
    @files = get_bams($dir);
	
	#parse results and determine if the pipeline should die, or warnings thrown
	my @results;
    @results = verify_percent_coverage( \@files, $size, $min_depth );
    
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
    my ( $files, $size, $min_depth ) = @_;

    # Get the number of available cores
    #----------------------------------#
    my $cores = get_num_cores();

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
        my $name = fileparse( $file, ".bam" );

        # Run samtools depth and get results
        #------------------------------------#
        my $result = `samtools depth $file`;

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

#-------------------------------------------#
# Get the number of cores - leave one alone #
#-------------------------------------------#
sub get_num_cores {
    my $num_cpus = `cat /proc/cpuinfo | grep processor | wc -l`;
    chomp $num_cpus;
    return $num_cpus - 1;
}

#-------------------------------------------#
# Create a File::Iterator for the bam files #
#-------------------------------------------#
sub open_dir {
    my ($dir) = @_;
    return new File::Iterator(
        DIR     => $dir,
        RECURSE => 0,
        FILTER  => sub { $_[0] =~ /\.bam$/ }
    );
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

#-----------------------------------------#
# Get array of all files in the directory #
#-----------------------------------------#
sub get_bams {
    my ($dir) =@_;
    
    opendir my ($dh), $dir;
    my @dirs = readdir $dh;
    closedir $dh;

    $dir =~ s/\/$//;

    #removing both '.' and '..' files and putting back the full path to a list
    return  sort { $a cmp $b  } grep { /\.bam$/ } map { "$dir/$_"} grep { !/^\.\.?/} @dirs;
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