#!/usr/bin/env perl
# consolidate_vcfs
# Purpose:  Given a pair of *.vcf files, consolidates them into one file.

use warnings;
use strict;
use FindBin;
use lib $FindBin::Bin.'/../lib';
use lib $FindBin::Bin;
use Pod::Usage;
use Getopt::Long;
# VCF module from vcftools: http://vcftools.sourceforge.net/index.html
use Vcf;
use File::Temp qw /tempdir/;
use File::Path qw /rmtree /;
use File::Copy;
use File::Basename;

__PACKAGE__->run unless caller;

my $verbose;


sub run
{
	my ($freebayes, $mpileup, $coverage_cutoff, $min_mean_mapping, $ao, $bcftools, $output, $filtered_density_out, $skip_density_filter, $window_size, $density_threshold) = prepare_inputs(@_);

	my $resulting_file = combine_vcfs($freebayes, $mpileup, $coverage_cutoff, $min_mean_mapping, $ao, $bcftools, $output, $filtered_density_out, $skip_density_filter, $window_size, $density_threshold);

	my $result = move($resulting_file,$output);

	if (not $result)
	{
		die "$!"
	}

	### Return values here, script ends

	return;
}


sub combine_vcfs{
    my ($freebayes, $mpileup, $coverage_cutoff, $min_mean_mapping, $ao, $bcftools, $output, $filtered_density_out, $skip_density_filter, $window_size, $density_threshold) = @_;
    
    #create temp working directory for combining the VCFs
    my $template = "consolidate_vcfs-XXXXXX";
    my $tmp_dir = tempdir ($template, TMPDIR=> 1, CLEANUP => 1);

    my ($ext,$out_type) = ('.bcf','-O b'); # intermediate files will be in binary format for speed ups (bcf.gz)
    my $cmd;

    my $file_name = "$tmp_dir/strain.bcf.gz";

    my ($dir) = "$tmp_dir/isec_dir";

    if ( not -d $dir)
    {
        mkdir $dir;
    }

    #before running anything else
    #need to run filtered-coverage on original mpileup bcf file
    #going to use the default one provided from bcftools filter instead of a custom one.
    #if we are using vcf files, need to convert first to vcf before applying the filter. The reason is bug with bcftools where SOME files will put the wrong FLAG in....
    if ( $ext eq '.vcf.gz')
    {
        $cmd = "$bcftools view $mpileup -O v | $bcftools  filter -s 'coverage' -i 'DP>=$coverage_cutoff || DP4[0]+DP4[1]>=$coverage_cutoff'  $out_type > $dir/coverage_mpileup$ext";
    }
    else {
        $cmd = "$bcftools  filter -s 'coverage' -i 'DP>=$coverage_cutoff || DP4[0]+DP4[1]>=$coverage_cutoff' $mpileup $out_type > $dir/coverage_mpileup$ext";
    }

    system($cmd) == 0 or die "Could not run $cmd";


    $cmd = "$bcftools index  $dir/coverage_mpileup$ext";
    system($cmd) == 0 or die "Could not run $cmd";


    #confirm SNVS in freebayes by comparing them to mpileup REF and ALT columns
    #mark all SNVS found in mpileup but NOT in freebayes as filtered-mpileup with 'some' option. 'some' options allows
    #only records where some subset of ALT alleles match are compatible
    #                  #chrom      #pos     #ref  #alt
    #so if mpileup had NC_007530.2|668709 . T     G,A
    #it will map to a freebayes with
    #                  NC_007530.2|668709 . T     G
    $cmd = "$bcftools  isec $freebayes $dir/coverage_mpileup$ext -p $dir -c some $out_type";
    system($cmd) == 0 or die "Could not run $cmd";

	#result from bcftools isec is a directory that contains multiple *$ext files
	#ignoring: 0000$ext is for unique position just for vcf-split aka freebayes (we do not use at all since should be empty because mpileup report ALL positions)
	#use: 0001$ext is for unique position only found in mpileup
	#use: 0002$ext is positions that both freebayes and mpileup have a consensus on the base pair call (either a SNV or same as reference)
	#VCF line that is kept is the one from freebayes and NOT mpileup
	#use: 0003$ext same as 0002$ext but where mpileup VCF line is kept and not freebayes. Need so we can confirm isec SNVS from freebayes (0002$ext)


    #filter out SNVs that were only found in mpileup but NOT in freebayes
    $cmd = "$bcftools  filter  -m + -s 'mpileup' -i ' TYPE!=\"snp\" ' $dir/0001$ext $out_type  > $dir/1-0001$ext";

    system($cmd) == 0 or die "Could not run $cmd";

	#need to get rid of the stupid format information since they cannot be merge later downstream
	$cmd ="$bcftools  view -h $dir/1-0001$ext";
	my $result = `$cmd`;
	if ($result =~ /##FORMAT=\<ID=GL/)
	{
	    $cmd = "$bcftools  annotate -x FORMAT -x FORMAT/GT -x FORMAT/GL  $dir/1-0001$ext $out_type > $dir/filtered_mpileup$ext";
	}
	elsif ( $result eq '')
	{
	    die "Failed to retrieve header of '1-0001$ext' for strain '$mpileup'\n";
	}
	else
	{
		$cmd = "$bcftools  annotate -x FORMAT -x FORMAT/GT $dir/1-0001$ext $out_type > $dir/filtered_mpileup$ext";
	}
	system($cmd) == 0 or die "Could not run $cmd";

	######################################################################################################


    #Doing two level of filtering here
    #first is by alternative allele ratio (snv abundance ratio). Default is 75% of the reads need to show a SNV exist
    #Second is MQM (min mean mapping) needs to be above a threshold (default is 30)
    #if either of them fail, it will be hard clip out.
    #NB that not sure how it handles when have multiple different alternative alleles
    #also hard clipping ones that fail filtering. Do not want to have them appear in the snvalign-positions since they never passed
    $cmd = "$bcftools  filter  -m + -e  'MQM<$min_mean_mapping || INFO/AO/INFO/DP<$ao'  $dir/0002$ext $out_type   > $dir/1-0002$ext && bcftools index $dir/1-0002$ext";
    system($cmd) == 0 or die "Could not run $cmd";



    my $mpileup_checked_bcf = check_reference($bcftools,"$dir/1-0002$ext","$dir/0003$ext",$dir,"$dir/filtered_freebayes$ext",$out_type,$ext);

    if ($mpileup_checked_bcf )
    {
        die "Could not correctly format intersection mpileup file\n";
    }


    $cmd = "$bcftools  annotate -x FORMAT -x FORMAT/GL -x FORMAT/GQ  $dir/filtered_freebayes$ext $out_type > $dir/filtered_freebayes2$ext";
    system($cmd) == 0 or die "Could not run $cmd";

    $cmd = "$bcftools index  $dir/filtered_freebayes2$ext";
    system($cmd) == 0 or die "Could not run $cmd";

    $cmd = "$bcftools index  $dir/filtered_mpileup$ext";
    system($cmd) == 0 or die "Could not run $cmd";

    #need to have mpileup header otherwise bcftools index has issues
    #reason being is that freebayes does not report ##contig which is needed for bcftools index/query. Issue arise with some edge cases
    $cmd = "$bcftools merge --print-header  $dir/filtered_mpileup$ext $dir/filtered_freebayes2$ext > $dir/header";
    system($cmd) == 0 or die "Could not run $cmd";

    $cmd = "$bcftools  merge -O b --use-header $dir/header $dir/filtered_freebayes2$ext $dir/filtered_mpileup$ext > $file_name";
    system($cmd) == 0 or die "Could not run $cmd";

    $cmd = "$bcftools index -f $file_name";
    system($cmd) == 0 or die "Could not run $cmd";

    if(not defined $skip_density_filter)
    {
        $cmd = "$bcftools plugin filter_snv_density $dir/filtered_freebayes2$ext  --  --region_file $filtered_density_out";
            
        
        if(defined $window_size)
        {
            $cmd .= " --window_size $window_size";
        }
            
        if(defined $density_threshold)
        {
            $cmd .= " --threshold $density_threshold";
        }
            
        system($cmd) == 0 or die "Could not run $cmd";
    }
    
    rmtree $dir;

    print STDERR "$1" if ($1);

    return $file_name;

}


sub check_reference
{
    my ($bcftools,$freebayes,$mpileup,$basedir,$output,$out_type,$ext) = @_;
    my $cmd;

    #check to see if we have any records to run again
    my $stats = `$bcftools  stats  $freebayes`;
    if ( $stats =~ /number of records:\s+(\d+)/)
    {
        #if we do not find any valid SNVs in freebayes files,
        #then we just use the original 0003 vcf file produced by bcftools isec
        my $filter_vcf = $mpileup;
        
        if ($1)
        {
            #check to ensure that there is no reference in the header that does NOT at least one have record.
            #If there is no record for a reference
            #bcftools will simply freeze up
            
            #get header line and parse out the ##contig=
            my $out = `$bcftools view -h $mpileup | grep "##contig" | sed -e 's/^##contig=<ID=//' -e 's/,length=[0-9][0-9]*>//'`;
            die "Error: no ##contig entries in '$out'" if ($out =~ //); #assume that we have at least one reference

            my %refs;

            foreach ( split /\n/,$out)
            {
                $refs{$_}++; #adding each reference found into a hash
            }

            #determine which references ACTUALLY have position(s) in the body of the vcf.
            $out = `$bcftools index -s  $mpileup`;
            die "Error: no entry found using bcftools index  in '$out'" if ($out =~ //); #assume that we have at least one reference coming back from bcftools index
            foreach my $line( split/\n/,$out)
            {
                my @data=split/\t/,$line;
                if (exists $refs{$data[0]}) #removing reference(s) that have position in the body
                {
                    delete $refs{$data[0]};
                }
            }

            #if any reference still in %refs, means there is NO position in the vcf and need to be removed!
            if ( keys %refs)
            {
                $filter_vcf = "$basedir/corrected_mpileup.vcf.gz"; #will need a new modified VCF file
                
                my $old_header = `$bcftools view -h $mpileup`; #grabbing the current header that needs to be cleaned of references that do NOT have positions
                open my $out, ">$basedir/newheader";

                foreach my $line( split/\n/,$old_header)
                {
                    my $filter_in=1;
                    foreach my $ref( keys %refs)
                    {
                        if ( (index $line,$ref) !=-1)
                        {
                            #line was found and will not be added to the new header file
                            $filter_in=0;
                            last; #no point looking farther since only one entry per header
                        }
                    }
                    #if true, then print out header line into new file
                    if ( $filter_in)
                    {
                        print $out "$line\n";
                    }
                }
                close $out;

                $cmd = "$bcftools view -H  $mpileup > $basedir/no_header_vcf";
                system($cmd) == 0 or die "Could not run $cmd";

                #combine new vcf header and the original body together as a vcf
                $cmd = "cat $basedir/newheader $basedir/no_header_vcf > $basedir/new_mpileup.vcf";
                system($cmd) == 0 or die "Could not run $cmd";

                #need to convert into compressed and index due to edge case where 0003
                #from bcftools isec mpileup may have extra contig entries which will
                #cause annotate command below to die
                $cmd = "$bcftools view $basedir/new_mpileup.vcf -O z > $filter_vcf";
                system($cmd) == 0 or die "Could not run $cmd";
                
                $cmd = "$bcftools index -f $filter_vcf";
                system($cmd) == 0 or die "Could not run $cmd";
            }

            #applying mpileup filter-coverage status to intersect SNVs of both freebayes and mpileup
            $cmd = "$bcftools  annotate  $freebayes -a $filter_vcf $out_type -c FILTER > $output";
            system($cmd) == 0 or die "Could not run $cmd";
        }
        else
        {
            $cmd = "ln -s $freebayes  $output";
            system($cmd) == 0 or die "Could not run $cmd";
        }
    }
    else
    {
        die "Could not run bcftools stats on file '$freebayes'\n";
    }

    return 0;
}


sub prepare_inputs
{
		my ( $man, $help, $freebayes, $mpileup, $coverage_cutoff, $min_mean_mapping, $ao, $bcftools, $output, $filtered_density_out, $skip_density_filter, $window_size, $density_threshold);

		if( @_ && $_[0] eq __PACKAGE__ )
		{
			GetOptions(
				"vcfsplit=s"               => \$freebayes,
				"mpileup=s"                => \$mpileup,
				"coverage-cutoff|c=i"      => \$coverage_cutoff,
				"min-mean-mapping=i"       => \$min_mean_mapping,
				"snv-abundance-ratio=s"    => \$ao,
				"b|bcftools-path=s"        => \$bcftools,
				"o|output=s"               => \$output,
				"f|filtered-density-out=s" => \$filtered_density_out,
				"s|skip-density-filter"    => \$skip_density_filter,
				"w|window-size=i"          => \$window_size,
				"d|density-threshold=i"    => \$density_threshold,
				"h|help"                   => \$help,
				"m|an"                     => \$man,
				"v|verbose"                => \$verbose
			);
			pod2usage(0) if $help;
			pod2usage(-verbose => 2, -exitval=> 0) if $man;
		}
		else
		{
			($freebayes, $mpileup, $coverage_cutoff, $min_mean_mapping, $ao, $bcftools, $output, $filtered_density_out, $skip_density_filter, $window_size, $density_threshold) = @_;
		}

		$verbose = 0 if (not defined $verbose);

		if(not defined $freebayes)
		{
			print STDERR "No freebayes file specified!\n\n";
			pod2usage(1);
		}
		if(not -e $freebayes)
		{
			print STDERR "Unable to find freebayes file '$freebayes'\n\n";
			pod2usage(1);
		}

        if (not defined $mpileup)
	    {
            print STDERR "No mpileup file specified!\n\n";
		    pod2usage(1);
        }
	    if (not -e $mpileup)
	    {
		    print STDERR "Unable to find mpileup file '$mpileup'\n\n";
		    pod2usage(1);
	    }

		if (not defined $output)
		{
			print STDERR "No output specified.";
			pod2usage(1);
		}

        if (not defined $bcftools)
        {
            #check to see if bcftools is on the path
            #normally we always want to be passed the path to the tool but this is for Galaxy implementation
            my $alive=`bcftools 2>&1 1>/dev/null`;
            if ( $alive && $alive =~ /Program: bcftools/)
            {
                $bcftools="bcftools";
            }
        else
        {
            print STDERR "bcftools-path not defined and not found on the path.\n";
			pod2usage(1);
        }
    }

    if(not defined $min_mean_mapping)
    {
    	$min_mean_mapping = 30;
    }

    #need to have ratio/percentage in double format i.e 75% = 0.75
    if(not defined $ao)
    {
       $ao = 0.75;
       print "Value for '--snv-abundance-ratio' not defined, defaulting to 0.75\n";
    }
    elsif ( $ao > 1)
    {
        print "Assuming that '$ao' is given as percentage for alternative allele. Changing to a decimal\n";
        $ao = $ao/100;
    }

    #need check to see if bcftools was complied with htslib and also has the correct plugin installed
    my $usage_state = `$bcftools 2>&1 1>/dev/null`;
    if ( not $usage_state =~ /Version: .* \(using htslib/ )
    {
        print STDERR "bctools was not complied with htslib.\nPlease re-compile with htslib\nInstruction: http://samtools.github.io/bcftools/\n";
				pod2usage(1);
    }

    if (not defined $coverage_cutoff)
    {
        print STDERR "warning: coverage-cutoff not set, assuming it is 1\n";
        $coverage_cutoff = 1;
    }
    elsif ($coverage_cutoff !~ /^\d+$/)
    {
        print STDERR "coverage-cutoff=$coverage_cutoff is invalid\n";
				pod2usage(1);
    }

	if (defined $skip_density_filter)
	{
		print "Skipping density filter";
	}
	else
	{
		if (not defined $filtered_density_out)
		{
			print STDERR "Filtered density output not specified!\n";
			pod2usage(1);
		}
		if (not defined $window_size)
		{
			print "Window size not specified. Using default value...";
		}
		if (not defined $density_threshold)
		{
			print "Density threshold not specified. Using default value...";
		}
	}

    return ($freebayes, $mpileup, $coverage_cutoff, $min_mean_mapping, $ao, $bcftools, $output, $filtered_density_out, $skip_density_filter, $window_size, $density_threshold);
}

1;

=pod

=head1 NAME

consolidate_vcfs.pl

=head1 SYNOPSIS

consolidate_vcfs.pl --vcfsplit v1=files/dataset1.dat --mpileup v1=files/dataset2.dat --coverage-cutoff 15 --min-mean-mapping 30 --snv-abundance-ratio 0.75 --bcftools-path /opt/bcftools/bcftools

=head1 OPTIONS

=over

=item B<--vcfsplit> [REQUIRED]

Multiple list of key/value pair. Multiple .gz files can be input.  Example with 3 gz files: --vcfsplit 'name=/path/vcf1.gz' --vcfsplit 'name=/path/vcf2.gz' --vcfsplit 'name=/path/vcf3.gz'.

=item B<--mpileup> [REQUIRED]

Multiple list of key/value pair. Multiple .gz files can be input.  Example with 3 gz files: --mpileup 'name=/path/vcf1.gz' --mpileup 'name=/path/vcf2.gz' --mpileup 'name=/path/vcf3.gz'.

=item B<--coverage-cutoff> [REQUIRED]

The cutoff for coverage to include a reference base (default: 1).

=item B<--min-mean-mapping> [REQUIRED]

Mean mapping quality of observed alternate alleles.

=item B<--snv-abundance-ratio> [REQUIRED]

Ratio of number of reads that supporting a variant compare to total reads for a particular position. Percentage for values > 1 (e.g., 75 = 75%), decimal for values <= 1 (e.g., 0.75 = 75%).

=item B<--bcftools-path> [REQUIRED]

Full path to BCFTools.

=item B<-h>, B<--help>

Displays the help screen.

=back

=head1 DESCRIPTION

Consolidates a given set of *.vcf files for use by vcf2snv_alignment.

=cut
