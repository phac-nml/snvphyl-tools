#!/usr/bin/env perl
# vcf2pseudoalignmnet
# Purpose:  Given a set of *.vcf files, examines all called SNPs and generates a 'pseudoalignment'
#  of all the high quality SNPs/base calls (to be used for further phylogenetic analysis).

use warnings;
use strict;
use FindBin;
use lib $FindBin::Bin.'/lib';
use lib $FindBin::Bin;
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;
use Streaming;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::SimpleAlign;
# VCF module from vcftools: http://vcftools.sourceforge.net/index.html
use Vcf;
use File::Temp qw /tempdir/;
use InvalidPositions;
use List::MoreUtils qw/all any firstidx/;
use File::Path qw /rmtree /;



my $verbose;
my $unknown_base = 'N';

my $snp_info = {'removed' => {'insertions' => 0, 'deletions' => 0, 'multi' => 0, 'other' => 0},
		'total' => 0,
		'positions' => 0,
		'snps' => {'filtered-coverage' => 0, 'filtered-mpileup' => 0, 'kept' => 0}
		};

sub usage
{
	"Usage: $0 --vcf-dir [vcf dir] --mpileup-dir [mpileup dir] --output-base [base of output alignment file]\n".
	"Parameters:\n".
	"\t--vcf-dir: The directory containing the vcf files.\n".
        "\t--vcfsplit: Multiple list of key/value pair  'name=path/to/vcf.gz'".          
	"\t--mpileup-dir: Directory containing the vcf files produced by 'samtools mpileup \$file | bcftools view -cg' (*.vcf).\n".
        "\t--mpileup: Multiple list of key/value pair  'name=path/to/vcf.gz'".
	"\t-o|--output-base:  The output base name for the alignment file(s)\n".
	"Options:\n".
	"\t-u|--uniquify:  Make the seq names unique and print mapping file to real names (for phylip format limitations)\n".
	"\t-r|--reference:  The name of the reference to use in the alignment (default: reference)\n".
	"\t-f|--format:  The format to output the alignment to, one of the Bio::AlignIO supported formats (default: fasta)\n".
	"\t-c|--coverage-cutoff:  The cutoff for coverage to include a reference base (default: 1)\n".
	"\t--invalid-pos: A TSV file that contains a list of range(s) (one per line) of CHROM\\tSTART_POS\\tEND_POS\\n".
	"\t--verbose:  More information printed\n".
	"\t--keep-ambiguous:  Keep ambiguous characters in alignment output file\n".
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

sub variant_info_to_hash
{
	my ($info_string) = @_;

	my @info = split(/;/,$info_string);
	my %info_hash = ();
	foreach my $info_entry (@info)
	{
		if (defined $info_entry and $info_entry ne '')
		{
			my @parts = split(/=/, $info_entry);
			if (defined $parts[0] and $parts[0] ne '')
			{
				if (@parts > 1)
				{
					$info_hash{$parts[0]} = $parts[1];
				}
				else
				{
					$info_hash{$parts[0]} = undef;
				}
			}
		}
	}

	return \%info_hash;
}



############
### MAIN ###
############

# maps format name to format file extension
my %valid_formats = ('fasta' => 'fasta', 'phylip' => 'phy', 'clustalw' => 'cl');

my ($vcf_dir, $mpileup_dir, $output_base, @formats, $reference, $coverage_cutoff);
my ($uniquify, $help, $requested_cpus, $invalid, $keep_ambiguous,$fasta,$bcftools);
my (%vcf_files, %mpileup_files);

my $command_line = join(' ',@ARGV);

if (!GetOptions('vcf-dir|d=s' => \$vcf_dir,
                'vcfsplit=s' => \%vcf_files,
		'mpileup-dir|b=s' => \$mpileup_dir,
                'mpileup=s' => \%mpileup_files,
		'format|f=s' => \@formats,
		'output-base|o=s' => \$output_base,
		'reference|r=s' => \$reference,
                'fasta=s' => \$fasta,
		'coverage-cutoff|c=i' => \$coverage_cutoff,
		'uniquify|u' => \$uniquify,
		'invalid-pos=s' => \$invalid,
		'help|h' => \$help,
		'keep-ambiguous' => \$keep_ambiguous,
                'numcpus=i' => \$requested_cpus,
                'b|bcftools-path=s' => \$bcftools,     
		'verbose|v' => \$verbose))
{
	die "Invalid option\n".usage;
}

print usage and exit(0) if (defined $help);
$verbose = 0 if (not defined $verbose);

if ( $vcf_dir and $mpileup_dir)
{
    die "vcf-dir does not exist\n".usage if (not -e $vcf_dir);

    die "mpileup-dir does not exist\n".usage if (not -e $mpileup_dir);
}
elsif ( scalar keys %vcf_files == 0 or scalar keys %mpileup_files ==0)
{
    die "Was not able to find any vcf files from freebayes and/or mpileup.";
}

die "output-base undefined\n".usage if (not defined $output_base);

die "bcftools-path not defined\n".usage if (not defined $bcftools or not -e $bcftools);

if (not defined $reference)
{
	print STDERR "reference name not defined, calling it 'reference'\n";
	$reference = 'reference';
}



if (defined $invalid)
{
    if ( ! -e $invalid)
    {
	die "Was given an invalid position file but could not locate it '$invalid'\n";
    }
}
else
{
    print STDERR "invalid position file not defined, Will ignore step\n";
}


$requested_cpus = 1 if (not defined $requested_cpus);
$uniquify = 0 if (not defined $uniquify);

if (@formats <= 0){
	print STDERR "warning: format not defined, assuming fasta\n";
	@formats = ("fasta");
}
else
{
	for my $format (@formats)
	{
		die "unrecognized format '$format', must be one of '".join(' ', keys %valid_formats),"'\n" if (not defined $valid_formats{$format});
	}
}

if (not defined $coverage_cutoff)
{
	print STDERR "warning: coverage-cutoff not set, assuming it is 1\n";
	$coverage_cutoff = 1;
}
elsif ($coverage_cutoff !~ /^\d+$/)
{
	die "coverage-cutoff=$coverage_cutoff is invalid\n".usage;
}


my $dh;
# fill table vcf_files with entries like if only provided the vcf-dir
#  vcf1 => dir/vcf1.vcf.gz
#  vcf2 => dir/vcf2.vcf.gz
if ( $vcf_dir) {
    opendir($dh, $vcf_dir) or die "error opening directory $vcf_dir: $!";
    %vcf_files = map { /^(.*)\.vcf\.gz$/; $1 => "$vcf_dir/$_"} grep { /\.vcf\.gz$/ } readdir($dh);
    closedir($dh);
}

die "No *.vcf.gz files found in $vcf_dir.  Perhas you need to compress and index with 'tabix' tools\n".
"Example: bgzip file.vcf; tabix -p vcf file.vcf.gz" if (keys(%vcf_files) <= 0);

my $total_samples = (keys %vcf_files);





if ( $mpileup_dir)
{
    # create table of mpileup files corresponding to input freebayes/variant vcf files
    # assumes files are named with the same prefix
    # ex vcf-dir/file1.vcf.gz and mpileup-dir/file1.vcf.gz    
    my $mpileup_table = create_mpileup_table(\%vcf_files, $vcf_dir, $mpileup_dir);

    if (not defined $mpileup_table)
    {
        die "Error: vcf-dir contains unmatched files in mpileup-dir";
    }
    
    else
    {
        %mpileup_files = %{$mpileup_table};
    }

}
else
{
    if (scalar keys %mpileup_files != scalar keys %vcf_files )
    {
        die "Error: vcfsplit contains uneven number compare to mpileup-files";
    }
    my $m_name= join ('',sort {$a cmp $b } keys %mpileup_files);
    my $v_name= join ('',sort {$a cmp $b } keys %vcf_files);
    if ( $m_name ne $v_name) {
        die "Error: vcfsplit contains unmatched files to mpileup-files";
    }
    
}

my $refs_info = refs_info($fasta);


my $invalid_pos;

if ($invalid)
{
	my $invalid_positions_parser = InvalidPositions->new;
	$invalid_pos = $invalid_positions_parser->read_invalid_positions($invalid);
}


#create temp working directory for all combines vcf files
#in future make them stay around...
my $tmp_dir = tempdir (CLEANUP => 1);
    
#combine the mpileup and freebayes vcf files together
my $files = combine_vcfs(\%vcf_files,\%mpileup_files, $coverage_cutoff,$invalid_pos,$bcftools,$tmp_dir,$requested_cpus);



my $valid_positions = "$output_base/pseudoalign-positions.tsv";
my $ya = filter_positions($files,$refs_info,$invalid_pos,$valid_positions,$requested_cpus);





my @samples_list = sort {$a cmp $b } keys %vcf_files;
# my %chromosome_align;
# my $unique_count = 1;
# my %name_map; # used to map sample name to other information
# my %sample_map; # keeps track of which samples have which unique ids (so we can properly increment unique_count)
my %total_positions_map; # keep track of total positions, and if valid/not
# for my $chromosome (keys %$vcf_data)
# {
# 	my ($alignment,$total_positions) = variants_alignment($vcf_data->{$chromosome}, $chromosome, $reference, \@samples_list, $mpileup_data, $coverage_cutoff,$invalid_pos,$keep_ambiguous);
# 	$total_positions_map{$chromosome} = $total_positions;
# 	for my $sample (sort {$a cmp $b} keys %$alignment)
# 	{
# 		next if (@{$alignment->{$sample}->{'positions'}} <= 0); # no alignments

# 		if (not exists $sample_map{$sample})
# 		{
# 			$sample_map{$sample} = "sample$unique_count";
# 			$unique_count++;
# 		}

# 		my $sample_id = $sample_map{$sample};
# 		my $sample_name = ($uniquify) ? $sample_id : $sample;

# 		if (not defined $chromosome_align{$sample})
# 		{
# 			$name_map{$sample_name} = "$sample\t$chromosome:".join('|',@{$alignment->{$sample}->{'positions'}});
			
# 			$chromosome_align{$sample} = {'header' => $sample_name};
# 			$chromosome_align{$sample}->{'data'} = $alignment->{$sample}->{'alignment'};
# 		}
# 		else
# 		{
# 			$name_map{$sample_name} .= " $chromosome:".join('|',@{$alignment->{$sample}->{'positions'}});

# 			$chromosome_align{$sample}->{'data'} .= $alignment->{$sample}->{'alignment'};
# 		}
# 	}
# }
# # print alignment
# my $aln = Bio::SimpleAlign->new();
# for my $sample (sort {$a cmp $b} keys %chromosome_align)
# {
# 	my $id = $chromosome_align{$sample}->{'header'};
# 	my $data = $chromosome_align{$sample}->{'data'};

# 	if ($keep_ambiguous)
# 	{
# 		die "error: SNP alignment for $sample contains an invalid character"
# 			if ($data =~ /[^ATCGN]/);
# 	}
# 	else
# 	{
# 		die "error: SNP alignment for $sample contains an invalid character"
# 			if ($data =~ /[^ATCG]/);
# 	}

# 	my $seq = Bio::LocatableSeq->new(-seq => $data, -id => $id);
# 	$aln->add_seq($seq);
# }

# # sets displayname for each sequence
# for my $seq_id ($aln->each_seq)
# {
# 	my $start = $seq_id->start;
# 	my $end = $seq_id->end;
# 	my $id = $seq_id->id."/$start-$end";
# 	$aln->displayname($id, $seq_id->id);
# }

# # check if alignment is flush
# die "Alignment blocks are not all of the same length" if (not $aln->is_flush());

# for my $format (@formats)
# {
# 	my $output_file = "$output_base.".$valid_formats{$format};
# 	my $io = Bio::AlignIO->new(-file => ">$output_file", -format => $format,-idlength=>30);
# 	$io->write_aln($aln);
# 	print STDERR "Alignment written to $output_file\n";
# }


# print snp stats
# print "# Command Line\n";
# print "# $command_line\n";
# print "# SNP statistics\n";
# print "# Processed $total_samples samples\n";
# print "# Total variant called SNPs processed: ".$snp_info->{'total'},"\n";
# print "#\tRemoved ".$snp_info->{'removed'}->{'insertions'}," insertions\n";
# print "#\tRemoved ".$snp_info->{'removed'}->{'deletions'}," deletions\n";
# print "#\tRemoved ".$snp_info->{'removed'}->{'multi'}," multi\n";
# print "#\tRemoved ".$snp_info->{'removed'}->{'other'}," other\n";
# print "# Total Valid Positions: ".$snp_info->{'positions'},"\n";
# print "# Total SNPs to process: ".($snp_info->{'positions'}*$total_samples)."\n";
# print "#\tSNPs called as N's:\n";
# print "#\t\tLow Coverage: ".$snp_info->{'snps'}->{'filtered-coverage'}."\n";
# print "#\t\tVariant/mpileup differences: ".$snp_info->{'snps'}->{'filtered-mpileup'},"\n";
# print "#\t\tInvalid Position based on user provided file: ".$snp_info->{'snps'}->{'filtered-invalid'},"\n" if $invalid_pos;
# print "#\tValid SNPs for analysis: ".$snp_info->{'snps'}->{'kept'},"\n";
# print "# Positions file in $valid_positions\n";

# # print other information
# print "#\n#AlnName\tSampleName\tPositions\n";
# for my $name (sort keys %name_map)
# {
# 	print "$name\t",$name_map{$name},"\n";
# }


    

# open(my $vfh, ">$valid_positions") or die "Could not open $valid_positions: $!";
# print $vfh "#Chromosome\tPosition\tStatus\tReference\t";
# my @samples_sorted_list = sort {$a cmp $b} @samples_list;
# print $vfh join("\t",@samples_sorted_list);
# print $vfh "\n";
# for my $chr (keys %total_positions_map)
# {
# 	my $positions = $total_positions_map{$chr};
# 	for my $pos (sort {$a <=> $b} keys %$positions)
# 	{
# 		my $samples = $positions->{$pos}->{'samples'};
# 		my $ref = $positions->{$pos}->{'ref'};
# 		print $vfh "$chr\t$pos\t".$positions->{$pos}->{'status'}."\t$ref\t";
# 		my $first = 1;
# 		die "error in total_positions_map, for $chr:$pos, not enough sample entries" if (@samples_sorted_list != scalar(keys %$samples));
# 		my $id = 0;
# 		for my $sample (sort {$a cmp $b } keys %$samples)
# 		{
# 			die "error: sample name $sample different from ".$samples_sorted_list[$id] if ($samples_sorted_list[$id] ne $sample);
# 			if ($first)
# 			{
# 				$first = 0;
# 				print $vfh $samples->{$sample};
# 			}
# 			else
# 			{
# 				print $vfh "\t".$samples->{$sample};
# 			}
# 			$id++;
# 		}
# 		print $vfh "\n";
# 	}
# }
# close($vfh);



exit;





sub combine_vcfs{
    my ($vcf_files,$mpileup_files, $coverage_cutoff,$invalid_pos,$bcftools,$tmp_dir,$cpus) = @_;

    my %files;

    my $pm;
    my $num_cpus=`cat /proc/cpuinfo | grep processor | wc -l`;
    chomp $num_cpus;
    #ensure that you user cannot request more threads then CPU on the machine
    if ( $cpus > $num_cpus) {
        $cpus = $num_cpus;
    }

    $pm=Parallel::ForkManager->new($cpus);
    
    $pm -> run_on_finish ( # called BEFORE the first call to start()
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $child_data) = @_;
            # retrieve data structure from child
            if (defined($child_data)) {  # children are forced to send anything
                my ($name) = keys %$child_data;
                $files{$name} = $child_data->{$name};
            } else {
                die "One or more vcf file did not produce any data!\n";
            }
        }
    );

    foreach my $sample( keys %$mpileup_files) {
        
        my $pid = $pm->start and next;
        my $f_file = $vcf_files->{$sample};
        my $m_file = $mpileup_files->{$sample};

        my ($cmd,$status);
        
        
        my $file_name = "$tmp_dir/$sample" . "_combined.vcf";

        my ($dir) = "$tmp_dir/$sample" . '_answer';

        #confirm SNPS in freebayes by comparing them to mpileup REF and ALT columns
        #mark all SNPS found in mpileup but NOT in freebayes as filtered-mpileup with 'some' option. 'some' options allows
        #only records where some subset of ALT alleles match are compatible
        #                  #chrom      #pos     #ref  #alt 
        #so if mpileup had NC_007530.2|668709 . T     G,A
        #it will map to a freebayes with
        #                  NC_007530.2|668709 . T     G
        #$cmd = "$bcftools  isec $f_file $ready_mpileup -p $dir -c some -O z";
        $cmd = "$bcftools  isec $f_file $m_file -p $dir -c some -O z";
        $status = system($cmd);


        #filter with C complied nml specific filtering, also removing all information in FORMAT column, otherwise we cannot merge farther down
        $cmd = "$bcftools  annotate -x FORMAT -p filter_mpileup:dp=$coverage_cutoff $dir/0001.vcf.gz -O z  > $dir/filtered_mpileup.vcf.gz";
        $status = system($cmd);

        
        #filter by coverage and ratio of 75% with alternative allele
        #also filter by MQM flag = minumum mean mapping quality with > 30
        #NB that not sure how it handles when have multiple different alternative alleles
        #also hard clipping ones that fail filtering. Do not want to have them appear in the pseudo-positions since they never passed
        $cmd = "$bcftools  annotate -x FORMAT -p filter_freebayes:dp=$coverage_cutoff:mqm=30:ao=75  $dir/0002.vcf.gz -O z  > $dir/filtered_freebayes.vcf.gz";
        $status = system($cmd);
        
        
        
        #combine header but ignore GL tag for freebayes
        # $cmd = "zgrep '#' $dir/finish_freebayes.vcf.gz | grep -e '=GL' -e 'CHROM' -v | tail -n+5 > $dir/header_1";
        # $status = system($cmd);

        # $cmd = "zgrep '#' $dir/finish_mpileup.vcf.gz> $dir/header_2";
        # $status = system($cmd);

        # $cmd = "cat $dir/header_1 $dir/header_2 > $dir/header";
        # $status = system($cmd);


        $cmd = "$bcftools index -f $dir/filtered_freebayes.vcf.gz";
        $status = system($cmd);


        $cmd = "$bcftools index -f $dir/filtered_mpileup.vcf.gz";
        $status = system($cmd);

        #get the fake merge vcf header
        my $header = $FindBin::Bin . '/fake_vcf_header/header';
        
        
        $cmd = "$bcftools  merge --use-header $header $dir/filtered_freebayes.vcf.gz $dir/filtered_mpileup.vcf.gz > $file_name 2>/dev/null";

        $status = system($cmd);
        
        #bgzip up and then index
        `bgzip $file_name`;
        my $bgzip = $file_name . ".gz";
        
        `$bcftools index -t -f $bgzip`;

        rmtree $dir;
        $pm->finish(0,{"$sample" =>$bgzip});        

    
    
    }
    
    $pm->wait_all_children;

    return \%files;
}


sub filter_positions {
    my ($files,$refs,$invalid_pos,$valid_positions,$cpus) = @_;

    my %results;
    my $pm;
    my $num_cpus=`cat /proc/cpuinfo | grep processor | wc -l`;
    chomp $num_cpus;
    #ensure that you user cannot request more threads then CPU on the machine
    if ( $cpus > $num_cpus) {
        $cpus = $num_cpus;
    }
    

    $pm=Parallel::ForkManager->new($cpus);
    
    $pm->run_on_finish ( # called BEFORE the first call to start()
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $child_data) = @_;
            # retrieve data structure from child
            if (defined($child_data)) {  # children are forced to send anything
                my ($name) = keys %$child_data;
                $results{$name} = $child_data->{$name};
            } else {
                die "One or more vcf file did not produce any data!\n";
            }
        }
    );
    my $parallel=1;
    
    my $bit_size = 100000;
    my $job_id=0;


    #print header file
    open my $out, '>', $job_id;
    print $out "#Chromosome\tPosition\tStatus\tReference\t";
    my @samples_list = sort {$a cmp $b } keys %$files;
    print $out join("\t",@samples_list);
    print $out "\n";
    close $out;
    
    
    foreach my $chrom( keys %$refs) {
        
        my ($start,$stop)=(0,0);
        my ($range);
        my $length = $refs->{$chrom}{'length'};
        
        if ( !$length) {
            die "Could not determine length of '$chrom' from provided reference fasta file\n";
        }
        
        
        while ($stop < $length ) {
            #inclusive range, we do not care if we go over b/c vcf query will only return what it can.
            $start = $stop +1;
            $stop = $stop+$bit_size;
            $range="$chrom:" . join('-',($start,$stop));

            $range = quotemeta $range;

            $job_id++;

            my $pid = $pm->start and next if $parallel;
            
            my $f_name = $job_id;
            open my $out, '>',$f_name;
            
            my $streamers = Streaming::create_streamers($files,$range,$job_id);
            my $cur_pos = $start;
            
            my @data = $streamers->($chrom,$cur_pos);


            #search for all entries to have "EOF" as their data.
            while (any { $_->{'status'} ne 'EOF' } @data) {

                #get current bp for current positions
                my $ref_bp = $refs->{$chrom}{'bps'}[$cur_pos];
                
                #if we do not have any SNP, we do not print at all!
                if ( any {  exists $_->{'alt'} && $_->{'alt'} ne '.' } @data ) {
                    my @line = ($chrom,$cur_pos);

                    #all columns are passed all cut off parameters
                    if ( all { $_->{'status'} eq 'PASS' }  @data) {

                        if ($invalid_pos && exists $invalid_pos->{"${chrom}_${cur_pos}"} ) {
                            push @line, 'filtered-invalid';
                        }
                        else {
                            push @line, 'valid';
                        }
                        
                        push @line,$ref_bp;
                        foreach ( @data) {
                            if ( $_->{'alt'} eq '.') {
                                push @line,$ref_bp;
                            }
                            else {
                                push @line,$_->{'alt'};
                                }
                        }
                        print $out join("\t",@line) . "\n";

                    }
                    #now we are working with position where there is at least one SNP and at least one filtered-*
                    elsif ( any { ( $_->{'status'} eq 'PASS' || $_->{'status'} eq 'filtered-invalid' )  && $_->{'alt'} ne '.'} @data ) {


                        if ($invalid_pos && exists $invalid_pos->{"${chrom}_${cur_pos}"} ) {
                            push @line, 'filtered-invalid';
                        }
                        else {
                            #always take filtered-invalid before anything else if it there.
                            my $index = firstidx {$_->{'status'} eq 'filtered-invalid'  } @data;

                            #if we did not find a filtered-invalid, find another none valid status
                            if ( $index == -1) {
                                $index = firstidx {$_->{'status'} ne 'PASS' } @data;
                            }
                            
                            if ( $data[$index]->{'status'} eq '' || $data[$index]->{'status'} eq 'EOF') {
                                push @line,'filtered-coverage';
                            }
                            else {
                                push @line,$data[$index]->{'status'};
                            }
                        }

                        
                        push @line,$ref_bp;
                        
                        foreach my $col ( @data) {
                            my $status = $col->{'status'};
                            
                            #if '', implies there is no reads covering that $cur_pos
                            if ( $status eq '' || $status eq 'EOF') {
                                push @line,'-';
                            }
                            else {
                                #if we have filtered-mpileup, we have inconsistent calles between variant callers
                                if ($status eq 'filtered-mpileup' ) {
                                    push @line,'N';
                                }
                                #have a position that passes the cut-off parameter. Either show the SNP or the reference
                                elsif ( $status eq 'PASS' || $status eq 'filtered-invalid') {
                                    if ( $col->{'alt'} eq '.') {
                                        push @line,$ref_bp;
                                    }
                                    else {
                                        push @line,$col->{'alt'};
                                    }
                                }
                                else {
                                    push @line,'-';
                                }
                                
                            }
                        }
                        
                        print $out join("\t",@line) . "\n";
                        
                        
                    }#end else
                    
                    
                }#end if
                
                
                
                $cur_pos++;
                @data = $streamers->($chrom,$cur_pos);
            }
            
            close $out;
            
            $pm->finish(0,{$job_id =>$f_name}) if $parallel;
            $results{$job_id}=$f_name if not $parallel;
        }
        
    
    }
    $pm->wait_all_children if $parallel;   

    #combine all results files in order
    my @cmd = "cat 0";
    foreach ( sort {$a <=> $b } keys %results) {
         push @cmd, $results{$_};
    }
    my $cmd = join(' ' , @cmd) . " > $valid_positions";
    
    `$cmd`;
    unlink "0";
    
    map { unlink $_ } values %results;
    
    

    return;
}

sub refs_info {
    my ($file) = @_;


    my %refs;
    my $in = Bio::SeqIO->new(-format=>'fasta',-file=>$file);
    while ( my $seq = $in->next_seq()) {
        $refs{$seq->display_id()}{'length'} = $seq->length;
        my @bps = ('X');
        
        map { push @bps, uc $_ } split //, $seq->seq;
        $refs{$seq->display_id()}{'bps'} = \@bps;
    }


    return \%refs;
    
}


