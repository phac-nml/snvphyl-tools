package Align::Nucmer;

use warnings;
use strict;

use File::Temp 'tempdir';
use Bio::SeqIO;
use File::Basename;
use Cwd qw(abs_path getcwd);

# Builds a new Align::Nucmer object
# Input
#	verbose: Turn on verbose information
sub new {
	my ($class,$verbose) = @_;

	my $self = {};
	bless($self,$class);

	$verbose = 0 if (not defined $verbose);
	$self->{'verbose'} = $verbose;

	return $self;
}

sub _fasta_lengths {
	my ($self,$file) = @_;

	my %lengths;
	my $in = Bio::SeqIO->new(-format=>'fasta',-file=>$file);

	while ( my $seq = $in->next_seq()) {
		$lengths{$seq->display_id} = $seq->length();
	}

	return \%lengths;
}

# aligns the given files with nucmer and parses the alignment
# input
#	reference:  The reference file to align
#	query: The query file to align
# output
#	bp:  A hash table of alignments
sub align_and_parse
{
	my ($self,$reference,$query) = @_;

	my $verbose = $self->{'verbose'};

	my %bp;
	my %alignments;

	die "error: reference file is not defind" if (not defined $reference);
	die "error: reference=$reference file does not exist" if (not -e $reference);
	die "error: query file is not defind" if (not defined $query);
	die "error: query=$query file does not exist" if (not -e $query);

	my $nucmer_dir = tempdir("ncumer.XXXXXX", TMPDIR => 1, CLEANUP => !$verbose);
	print STDERR "Nucmer dir: $nucmer_dir\n" if ($verbose);
	my $cwd = getcwd;

	my $reference_name = basename($reference);
	my $query_name = basename($query);

	my $reference_path = abs_path($reference);
	my $query_path = abs_path($query);

	chdir($nucmer_dir);
	my $delta_prefix = "${reference_name}_${query_name}";
	my $delta_file = "$nucmer_dir/$delta_prefix.delta";
	my $command = "nucmer --prefix=$delta_prefix $reference_path $query_path 2> /dev/null 1> /dev/null";
	system($command) == 0 or die "Could not execute '$command'";

	my $delta_filter = "$nucmer_dir/$delta_prefix.filter.delta";
	$command = "delta-filter -q -r $delta_file > $delta_filter";
	system($command) == 0 or die "Could not execute '$command'";

	my $snvs_file = "$nucmer_dir/snvs";
	$command = "show-snps -CTr $delta_filter > $snvs_file";
	system($command) == 0 or die "Could not execute '$command'";

	chdir($cwd);

	my $ref_lengths = $self->_fasta_lengths($reference_path);
	my $query_lengths = $self->_fasta_lengths($query_path);

	for my $ref_id (keys %$ref_lengths)
	{
		for my $query_id (keys %$query_lengths)
		{
			my $align_file = "$nucmer_dir/${ref_id}_${query_id}.align";
			$command = "show-aligns $delta_filter $ref_id $query_id 1> $align_file 2>/dev/null";

			# shift by 8 to get real return value
			my $ret_value = (system($command) >> 8);
			if ($ret_value != 0)
			{
				if ($ret_value == 1)
				{
					print STDERR "no alignments for $ref_id $query_id, command '$command'" if ($verbose);
				}
				else
				{
					die "Error: could not execute '$command', retvalue=$ret_value";
				}
			}

			if (not exists $bp{$ref_id})
			{
				$bp{$ref_id} = undef;
			}

			$self->_parse_alignments(\%bp,$align_file);
		}
	}

	# remove any empty datasets
	for my $ref_id (keys %bp)
	{
		my $positions = $bp{$ref_id};
		if (keys %$positions <= 0)
		{
			# delete empty positions hash
			$bp{$ref_id} = undef;
		}
	}

	$alignments{'all'} = \%bp;
	$alignments{'snvs'} = $self->_handle_show_snvs($snvs_file);

	return \%alignments;
}

# returns a data structure containing SNVs, or undef if no SNVs
sub _handle_show_snvs
{
	my ($self,$show_snvs_file) = @_;

	my %snvs;

	die "error: show_snvs_file is undefined" if (not defined $show_snvs_file);
	die "error: show_snvs_file does not exist" if (not -e $show_snvs_file);

	open(my $fh, "<$show_snvs_file") or die "Could not open $show_snvs_file: $!";
	
	my $line;
	# skip first line
	$line = readline($fh);

	return undef if (not defined $line);
	$line = readline($fh);
	die "error: $show_snvs_file not valid show-snvs file. try running 'show-snvs -CTr'"
		if ($line ne "NUCMER\n");

	# skip next two lines
	readline($fh);
	readline($fh);

	while($line = readline($fh))
	{
		chomp $line;
		my @fields = split(/\t/,$line);
		my $ref_pos = $fields[0];
		die "error parsing $show_snvs_file, ref_pos undefined" if (not defined $ref_pos);

		my $ref = $fields[1];
		die "error parsing $show_snvs_file, ref not defined" if (not defined $ref);

		my $alt = $fields[2];
		die "error parsing $show_snvs_file, alt not defined" if (not defined $alt);

		my $ref_contig = $fields[8];
		die "error parsing $show_snvs_file, ref_contig not defined" if (not defined $ref_contig);

		$snvs{$ref_contig}->{$ref_pos} = {'ref' => $ref, 'alt' => $alt};		
	}

	close($fh);

	if (keys %snvs <= 0)
	{
		return undef;
	}
	else
	{
		return \%snvs;
	}
}

# Parses alignments and returns a hash table storing all positions
# input
#	show-aligns-file: The show-aligns file
#	show-snvs-file: The show-snvs file
# Output
#	bp:  A hash table storing all positions
sub _parse_alignments {
	#grabbing arguments either from command line or from another module/script
	my ($self,$bp,$show_aligns_file) = @_;

	my $verbose = $self->{'verbose'};
	
	my %bad;
	
	open my $in , '<', $show_aligns_file;

	#get skip first 3 lines in the file to ease
	for ( 0..2) {
		<$in>;
	}
	my $alignment = <$in>;
	return undef if (not defined $alignment);

	chomp $alignment;
	my ($ref,$alt) = $alignment =~ /.*between (.+) and (.+)/;
		
	#skip blank line
	<$in>;
		
	my $next_align = $self->_next_alignment($in);
	while ( my ($details) = $next_align->()) {
		if (! %$details) {
			last;
		}

		#go thru each base pair that was alignmned to the reference and record the position
		my ($qseq,$hseq) = ($details->{'query_seq'},$details->{'hit_seq'});

		#check to see if we are increasing or decreasing
		my ($pos,$next) =($details->{'start_hit'},1);
			
		if ( $details->{'orient_hit'} eq '-1') {
			$next=-1;
			$pos = $details->{'stop_hit'};
		}
			
		my @qseq = split//,$qseq;
		foreach my $ref_bp( split //,$hseq) {
			# keep consistent case
			$ref_bp = uc($ref_bp);
			my $query_bp = uc(shift @qseq);
			#we do not care about indels at the moment.
			if ( $query_bp eq '.') {
				print STDERR "Found deletion on query @ ref:'$ref:$pos'. Skipping\n" if $verbose;
				$pos +=$next;
			}
			elsif ($ref_bp eq '.') {
				print STDERR "Found insertion on query @ ref:'$ref:$pos'. Skipping\n" if $verbose;
				# no increment of position as we want to keep same reference coordinates
			}
			else {
				if ( exists $bp->{$ref}->{$pos}) {
					print STDERR "Seen already '$ref:$pos'. Removing.\n" if $verbose;
					delete $bp->{$ref}->{$pos}; # get rid of position already in the good pile
					$bad{$ref}{$pos}++; # ensure that if we see that position again that we ignore it
					$pos +=$next;
				}
				elsif ( exists $bad{$ref}{$pos}) {
					$pos +=$next;
				}
				else {		
					$bp->{$ref}->{$pos}={'alt'=>$query_bp,'ref'=> $ref_bp};
					#increment/decrement to next pos
					$pos +=$next;
				}
			}
		}
	}
	
	close $in;
}


sub _next_alignment {
	my ($self,$in) = @_;
	my $lines;
	
	return sub {
		local $/ = "END ";	
		$lines= <$in>;
		my %details;
		{
			my ($align);
			if ($lines) {
				local $/ = "\n";
				open $align ,'<',\$lines;
				my $header = <$align>;
				chomp $header;
				
				#if we see 'END alignment'... we will skip the line since the alignment are not evenly split since they do NOT have unique end markers
				#having unique end markers would make parser a lot easier. i.e '//' in genbank files or '>' in fasta files
				if ( $header =~ /^alignment.*/) {
					$header = <$align>;
					chomp ($header);
				}
				
				if ($header) {
					#print "$header\n";
					#get orientation,start & stop for both query and reference
					my @line = split/\s+/,$header;
					$details{'start_query'} =$line[10];
					$details{'stop_query'} =$line[12];
					$details{'orient_query'} =$line[9];
					
					$details{'start_hit'} =$line[5];
					$details{'stop_hit'} =$line[7];
					$details{'orient_hit'} =$line[4];
					$details{'header'} = $header;
				}

				while (my $hit_line = <$align>) {
					chomp $hit_line;
					
					if ($hit_line && $hit_line =~ /\d+\s+.*/) {
						#remove the whitespace and coordinate from the beginning of the line and just save the base pair
						my (undef,$hit_seq) = split/\s+/,$hit_line;
						
						#hit line always follows a query line
						my $query_line = <$align>;
						chomp $query_line;
						my (undef,$query_seq) = split/\s+/,$query_line;
						
						my $match_line = <$align>;
						
						if (length $query_seq == length $hit_seq ) {
							$details{'query_seq'} .= $query_seq;
							$details{'hit_seq'} .= $hit_seq;
						}
						else {
							die "Alignments are not the same length for HSP : '$header'\n";
						}
					}
				}
			}
		}
		
		return \%details;
	}
}

1;
