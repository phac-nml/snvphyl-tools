package InvalidPositions;

use warnings;
use strict;

sub new
{
	my ($class) = @_;

	my $self = {};
	bless($self,$class);

	return $self;
}

# read the invalid positions file and returns a table of the positions
# input
#	invalid_positions_file: The file to read
# output
#	a table containing the invalid positions, formatted like
#		{ chrom_pos => 1 } if chromosome_position is invalid
#       a table that containing total number of invalid positions
#               { chrom => total_invalid }
sub read_invalid_positions
{
	my ($self,$invalid_positions_file) = @_;

	my %invalid;
        my %count_per_ref; #keep track of positions that were marked as invalid. 
        my %invalid_per_ref; #summary of %count_per_ref, we only want to return how many positions were marked as invalid. 
         

	open(my $fh, "<" , "$invalid_positions_file") or die "Could not open $invalid_positions_file: $!";

	while(my $line = readline($fh))
	{
		chomp $line;
		my ($sub_line) = ($line =~ /^([^#]*)/);
		next if ($sub_line eq ''); # comment line

		my ($chrom,$start,$end) = split(/\s+/,$sub_line);
		if (not defined $chrom or $chrom eq '')
		{
			print STDERR "warning: line '$line' contains no information for invalid positions file $invalid_positions_file, skipping...\n";
			next;
		}
		if ($start !~ /^\d+$/)
		{
			print STDERR "Warning: line '$line' not valid for invalid positions file $invalid_positions_file\n";
			next;
		}
		if ($end !~ /^\d+$/)
		{
			print STDERR "Warning: line '$line' not valid for invalid positions file $invalid_positions_file\n";
			next;
		}

		# swap in case start/end are reversed
		my $real_start = ($start < $end) ? $start : $end;
		my $real_end = ($start < $end) ? $end : $start;


		foreach my $i ( $real_start..$real_end ) {
			$invalid{"${chrom}_${i}"} = 1;
                        $count_per_ref{$chrom}{$i}++;
		}
	}

	close($fh);


        #summaries up how many positions were marked as invalid per reference
        foreach my $ref ( keys %count_per_ref) {
            $invalid_per_ref{$ref} = scalar keys %{$count_per_ref{$ref}};
        }
        
	return  (\%invalid,\%invalid_per_ref);
}

1;
