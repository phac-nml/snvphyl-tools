#!/usr/bin/env perl
# vcf2core
# Purpose:  
use warnings;
use strict;

use Getopt::Long;
use Storable qw /dclone store retrieve/;
use File::Basename;
use Parallel::ForkManager;
use FindBin qw($Bin);
use lib "$Bin";
# VCF module from vcftools: http://vcftools.sourceforge.net/index.html
use Vcf;
use File::Temp qw /tempdir tempfile/;
use List::MoreUtils qw/all/;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;

my $verbose;

sub usage
{
	"Usage: $0 --mpileup-dir [mpileup dir] --output-base [base of output alignment file]\n".
	"Parameters:\n".
	"\t--mpileup-dir: Directory containing the vcf files produced by 'samtools mpileup \$file | bcftools view -cg' (*.vcf).\n".
	"\t-o|--output-base:  The output base name for the alignment file(s)\n".
	"Options:\n".
	"\t-c|--coverage-cutoff:  The cutoff for coverage to include a reference base (default: 1)\n".
	"\t--verbose:  More information printed\n".
        "\t--gview_path: Full path to gview.jar".
        "\t--gview_style: Full path to default gview stylesheet".
	"\t-h|--help:  Help\n";
}

sub create_mpileup_table
{
	my ($vcf_files, $vcf_dir, $mpileup_dir) = @_;
	my %mpileup_table;

	for my $vcf_name (keys %$vcf_files)
	{
		my $vcf_file = $vcf_files->{$vcf_name};
		my $mpileup_file = "$mpileup_dir/".basename($vcf_file);
		$mpileup_table{$vcf_name} = $mpileup_file;
		if (not (-e $mpileup_file))
		{
			print STDERR "Could not find mpileup file \"$mpileup_file\" corresponding to vcf-file \"$vcf_file\"\n";
			return undef;
		}
	}

	return \%mpileup_table;
}



############
### MAIN ###
############



my ($mpileup_dir,$output_base,$coverage_cutoff,$help,$requested_cpus,$fasta,$gview,$gview_style);

my $command_line = join(' ',@ARGV);

if (!GetOptions(
		'mpileup-dir|b=s' => \$mpileup_dir,
		'output-base|o=s' => \$output_base,
		'coverage-cutoff|c=i' => \$coverage_cutoff,
                'i|fasta=s' => \$fasta,
		'help|h' => \$help,
                'numcpus=i' => \$requested_cpus,
                'gview_path=s' => \$gview,
                'gview_style=s' => \$gview_style,
		'verbose|v' => \$verbose))
{
	die "Invalid option\n".usage;
}

print usage and exit(0) if (defined $help);
$verbose = 0 if (not defined $verbose);

die "mpileup-dir undefined\n".usage if (not defined $mpileup_dir);
die "mpileup-dir does not exist\n".usage if (not -e $mpileup_dir);

die "output-base undefined\n".usage if (not defined $output_base);

die "gview undefined\n".usage if (not defined $gview);
die "gview_style undefined\n".usage if (not defined $gview_style);


$requested_cpus = 1 if (not defined $requested_cpus);


    
if (not defined $coverage_cutoff)
{
	print STDERR "warning: coverage-cutoff not set, assuming it is 1\n";
	$coverage_cutoff = 1;
}
elsif ($coverage_cutoff !~ /^\d+$/)
{
	die "coverage-cutoff=$coverage_cutoff is invalid\n".usage;
}


my %mpileup_files;

my $dh;
# fill table mpileup_files with entries like
#  vcf1 => dir/vcf1.vcf.gz
#  vcf2 => dir/vcf2.vcf.gz
opendir($dh, $mpileup_dir) or die "error opening directory $mpileup_dir: $!";
%mpileup_files = map { /^(.*)\.vcf\.gz$/; $1 => "$mpileup_dir/$_"} grep { /\.vcf\.gz$/ } readdir($dh);
closedir($dh);

die "No *.vcf.gz files found in $mpileup_dir.  Perhas you need to compress and index with 'tabix' tools\n".
"Example: bgzip file.vcf; tabix -p vcf file.vcf.gz" if (keys(%mpileup_files) <= 0);


my $info = determine_core($output_base,\%mpileup_files,$requested_cpus,$fasta,$coverage_cutoff);


printf("#Percentage of base pair in the core: %.2f \n",($info->{'core'}/$info->{'length'}*100));


create_figures($gview,$gview_style,$info->{'gffs'});

exit;


sub determine_core
{
    my ($output_base,$mpileup_files,$requested_cpus,$fasta,$coverage_cutoff) = @_;

    my $pm;
    my $num_cpus=`cat /proc/cpuinfo | grep processor | wc -l`;
    chomp $num_cpus;
    #ensure that you user cannot request more threads then CPU on the machine
    if ( $requested_cpus > $num_cpus) {
        $requested_cpus = $num_cpus;
    }

    my ($total_core,$total_length);
    my %files;
    my @gffs;

    $pm=Parallel::ForkManager->new($requested_cpus);
    # data structure retrieval and handling
    $pm -> run_on_finish ( # called BEFORE the first call to start()
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $child_data) = @_;
            # retrieve data structure from child
            if (defined($child_data)) {  # children are forced to send anything
                $total_core += $child_data->{'core'};
                if ( exists $files{$child_data->{'chrom'}}) {
                    push @{$files{$child_data->{'chrom'}}},$child_data->{'file'};
                }
                else {
                    $files{$child_data->{'chrom'}}=  [($child_data->{'file'})];
                }

            } else {
                die "One or more vcf file did not produce any data!\n";
            }
        }
    );


    my %ref_length = %{ get_lengths($fasta)};


    my $bit_size = 100000;
    foreach my $chrom( keys %ref_length) {
        
        my ($start,$stop)=(0,0);
        my ($range);
        my $length = $ref_length{$chrom};
        if ( !$length) {
            die "Could not determine length of '$chrom' from provided reference fasta file\n";
        }
        $total_length+=$length;
        
        
        while ($stop < $length ) {
            #inclusive range, we do not care if we go over b/c vcf query will only return what it can.
            $start = $stop +1;
            $stop = $stop+$bit_size;
            $range="$chrom:" . join('-',($start,$stop));
            $range = quotemeta $range;
            $pm->start and next;
            my ($fh,$filename) = tempfile();
            my $gffout =  Bio::Tools::GFF->new(-file => ">$filename" ,
                                               -gff_version => 3);
            
            my $core=0;
            my $streamers = create_streamers($mpileup_files,$range);
            my $cur_pos = $start;
        
        
            my @data = $streamers->($chrom,$cur_pos);
            #search for all entries to have "EOF" as their data.
            #not sure how it will handle multiple chromosome... just yet
            my @ranges;
            my $last_range_start;
            while (all { $_ ne 'EOF' } @data) {
            
            
                if (scalar @data >=1 && all { $_ && $_->{'cov'} >= $coverage_cutoff } @data ) {
                    $core++;
                    $last_range_start=$cur_pos if ! $last_range_start;
                }
                else {
                    #end of a range, may be one base pair or could be hundreds
                    
                    if ( $last_range_start) {
                        write_range($last_range_start,$cur_pos-1,$gffout);
                    }
                    $last_range_start=0;
                }
                
                $cur_pos++;
                @data = $streamers->($chrom,$cur_pos);
            }
            
            

            if ( $last_range_start) {
                write_range($last_range_start,$cur_pos-1,$gffout);
            }
        

            
            $pm->finish(0,{'core' =>$core,'file'=>$filename,'chrom'=>$chrom});
        }
    
        $pm->wait_all_children;

        foreach my $chrom( keys %files) {
            my ($fh,$combine) = tempfile();

            #get the reference name and make it at least usable as a filename
            my $name = $chrom;
            $name =~ s/\|$//;
            $name =~ s/\|/_/g;
            my $final = "$output_base/$name" . '.gff';


            #combine all segment for a single reference

            foreach ( @ {$files{$chrom}}) {
                `sed 1d $_ >> $combine`
            }
            #sort based on the position
            `sort -n -k 4 -t "\t" $combine > $final`;
            push @gffs,$final;
        }
       
    }

    
    my %info;
    $info{'core'} = $total_core;
    $info{'length'} = $total_length;
    $info{'gffs'} = \@gffs;
    
    
    return \%info;
}



sub create_figures
{
    my ($gview,$style,$gffs) = @_;
    
    
    #create gview image using all the gffs provided
    foreach my $gff( @$gffs) {
        my $out = $gff;
        $out =~ s/\.gff$/.png/;
        
        #time java -Xmx12G  -jar gview.jar -i task_1/NC_003997.3.gb -s style.gss -l circular -f png -o image.png -g results.gff -z 4
        `java -Xmx4G -jar $gview -i $fasta -s $style -l circular -f png -o $out -g $gff 15300 -H 4500   -z 4`;
    }
    
    return;
}
    

sub write_range {
    my ($x,$y,$gffout) = @_;
    #print "$x-$y\n";
    
    my $f = Bio::SeqFeature::Generic->new(
        -primary    => 'region',
        -seq_id     => 'task_2',
        -start      => $x,
        -end        => $y,
        -strand     => 1,
        -score      => 100.00,
        -frame      => 0,
    );
    $gffout->write_feature($f);
    
    return;
}

sub get_lengths
{
    my ($fasta) = @_;

    my %lengths;
    
    my $in = Bio::SeqIO->new(-format=>'fasta',-file=>$fasta);

    while ( my $seq = $in->next_seq())
    {
        $lengths{$seq->display_id} = $seq->length;
    }
    return \%lengths;
}



sub create_streamers {
    my ($files,$range) = @_;
    
    my @fhs;
    foreach ( sort {$a cmp $b } keys %$files ) {
        my %info;
        #fetch only the range from the reference that was assigned to this filehandle.
        my $cmd = "htscmd  vcfquery -r $range -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/DP\n' " . $files->{$_};
        
        #store results
        my $result = `$cmd`;
        die "'$cmd' produce no result\n" if ! $result;
        
        #created filehandle from scalar so we do not have to write to disk
        open (my $OUTPUT,'<',\$result);
        $info{'vcf'} = Vcf->new(fh => $OUTPUT ,_version_set=> '4.0');
        push @fhs,\%info;
    }
    
    return sub {
        my ($chrom,$pos) = @_;
        my %line;
        my @result;

        foreach my $info( @fhs) {
            while (1) {
                #check to see if we already readed a line but did not process it
                my $data = $info->{'last'} ? $info->{'last'} : $info->{'vcf'}->next_data_array;
                
                $info->{'last'} = undef if $info->{'last'};
                
                #if undef, we are at the end of the 'file'
                if (! $data){
                    push @result,'EOF';
                    last;
                }
                else{
                    my $chrom_cur = $data->[0];
                    my $position_cur = $data->[1];
                    
                    # found a position to keep
                    if ($chrom eq $chrom_cur && $pos == $position_cur){

                        my $ref = $data->[3];
                        my $alt = $data->[4];
                        my $coverage = $data->[7]; #atm just the DP in that column
                        push @result,{'ref' => $ref, 'alt' => $alt, 'cov' => $coverage};
                        last;
                    }
                    elsif ( $chrom eq $chrom && $pos < $position_cur) {
                        push @result,'';
                        $info->{'last'}= $data;
                        last;
                    }
                }
            }
                
        
            
        }
        
        return @result

    };
    
}

