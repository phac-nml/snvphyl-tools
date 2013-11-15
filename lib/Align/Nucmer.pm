package Align::Nucmer;

use warnings;
use strict;

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

# aligns the given files with nucmer and parses the alignment
# input
#	reference:  The reference file to align
#	query: The query file to align
# output
#	bp:  A hash table of alignments
#sub align_and_parse
#{
#	my ($self,$reference,$query) = @_;
#
#	die "error: reference file is not defind" if (not defined $reference);
#	die "error: reference=$reference file does not exist" if (not -e $reference);
#	die "error: query file is not defind" if (not defined $query);
#	die "error: query=$query file does not exist" if (not -e $query);
#
#	# do alignment here
#	my $show_aligns_file;
#	my $show_snps_file;
#
#	return $self->parse_alignments($show_aligns_file,$show_snps_file);
#}

sub _handle_show_snps
{
	my ($self,$show_snps_file) = @_;

	my %snps;

	die "error: show_snps_file is undefined" if (not defined $show_snps_file);
	die "error: show_snps_file does not exist" if (not -e $show_snps_file);

	open(my $fh, "<$show_snps_file") or die "Could not open $show_snps_file: $!";
	
	my $line;
	# skip first line
	$line = readline($fh);

	$line = readline($fh);
	die "error: $show_snps_file not valid show-snps file. try running 'show-snps -CTr'"
		if ($line ne "NUCMER\n");

	# skip next two lines
	readline($fh);
	readline($fh);

	while($line = readline($fh))
	{
		chomp $line;
		my @fields = split(/\t/,$line);
		my $ref_pos = $fields[0];
		die "error parsing $show_snps_file, ref_pos undefined" if (not defined $ref_pos);

		my $ref = $fields[1];
		die "error parsing $show_snps_file, ref not defined" if (not defined $ref);

		my $alt = $fields[2];
		die "error parsing $show_snps_file, alt not defined" if (not defined $alt);

		my $ref_contig = $fields[8];
		die "error parsing $show_snps_file, ref_contig not defined" if (not defined $ref_contig);

		$snps{$ref_contig}->{$ref_pos} = {'ref' => $ref, 'alt' => $alt};		
	}

	close($fh);

	return \%snps;
}

# Parses alignments and returns a hash table storing all positions
# input
#	show-aligns-file: The show-aligns file
#	show-snps-file: The show-snps file
# Output
#	bp:  A hash table storing all positions
sub parse_alignments {
	#grabbing arguments either from command line or from another module/script
	my ($self,$show_aligns_file,$show_snps_file) = @_;

	my $verbose = $self->{'verbose'};
	
	my %alignments;
	my %bp;
	my %bad;
	
	open my $in , '<', $show_aligns_file;

	#get skip first 3 lines in the file to ease
	for ( 0..2) {
		<$in>;
	}
	my $alignment = <$in>;
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
			if ( $query_bp eq '.' or $ref_bp eq '.') {
				print "Found indel @ '$pos'. Skipping\n" if $verbose;
				$pos +=$next;
				next;
			}
			else {
				if ( exists $bp{$ref}->{$pos} && $bp{$ref}->{$pos}->{'alt'} ne $query_bp) {
					print "Seen already '$pos' with " . $bp{$ref}->{$pos}->{'alt'} ." against $query_bp. Removing both entries\n" if $verbose;
					delete $bp{$ref}->{$pos}; # get rid of position already in the good pile
					$bad{$ref}{$pos}++; # ensure that if we see that position again that we ignore it
					$pos +=$next;
					next;
				}
				elsif (exists $bp{$ref}->{$pos} && $bp{$ref}->{$pos}->{'alt'} eq $query_bp ) {
					print "Same base pair for $pos\n" if $verbose;
				}
				elsif ( exists $bad{$ref}{$pos}) {
					$pos +=$next;
					next;
				}
				
				$bp{$ref}->{$pos}={'alt'=>$query_bp,'ref'=> $ref_bp};
				#increment/decrement to next pos
				$pos +=$next;
				
			}
			
		}
			
	}
	
	close $in;

	$alignments{'all'} = \%bp;
	$alignments{'snps'} = $self->_handle_show_snps($show_snps_file);
	
	return \%alignments;
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
