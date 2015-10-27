package CorePositions;

use warnings;
use strict;

use File::Basename;
use InvalidPositions;

# creates a new core positions parser
# input
#	core_file: the core positions file
#	bad_pos_file: the bad positions file
#	core_format: the core format (either 'gff' or 'tsv')
sub new
{
	my ($class, $core_file, $bad_pos_file, $core_format) = @_;

	my $self = {};
	bless($self,$class);

	die "error: core_format invalid" if ($core_format ne 'tsv' and $core_format ne 'gff');
	die "error: core_file is not defined" if (not defined $core_file);
	die "error: $core_file does not exist" if (not -e $core_file);
	if ($core_format eq 'gff')
	{
		die "error: $core_file not a proper gff file" if ($core_file !~ /\.gff$/i);
	}
	elsif ($core_format eq 'tsv')
	{
		die "error: $core_file not a proper tsv file" if ($core_file !~ /\.tsv$/i);
	}

	die "error: bad_pos_file is not defined" if (not defined $bad_pos_file);
	die "error: $bad_pos_file does not exist" if (not -e $bad_pos_file);

	my ($core_chrom,$core_pos) = $self->_parse_core_positions($core_file, $core_format);
	die "error: core chromosome is undefined" if (not defined $core_chrom);

	$self->{'core_pos_count'} = scalar(keys %$core_pos);

	$self->_remove_bad_pos($core_chrom,$core_pos,$bad_pos_file);

	$self->{'used_pos_count'} = scalar(keys %$core_pos);

	$self->{'core_file'} = $core_file;
	$self->{'bad_pos_file'} = $bad_pos_file;
	$self->{'core_pos_table'} = $core_pos;

	return $self;
}

sub get_core_pos_count
{
	my ($self) = @_;

	return $self->{'core_pos_count'};
}

sub get_used_pos_count
{
	my ($self) = @_;

	return $self->{'used_pos_count'};
}

sub get_core_positions
{
	my ($self) = @_;

	return $self->{'core_pos_table'};
}

sub _parse_coords
{
	my ($self,$chrom,$start,$end,$line) = @_;

	my ($real_start, $real_end) = (undef, undef);
	# swap in case start/end are reversed
        if (defined $chrom and $chrom ne '' and defined $start and defined $end
                and $start =~ /^\d+$/ and $end =~ /^\d+$/)
        {
                $real_start = ($start < $end) ? $start : $end;
                $real_end = ($start < $end) ? $end : $start;
        }
	else
	{
		die "error with start and end for line: $line";
	}

	return ($real_start, $real_end);
}

sub _parse_gff_line
{
	my ($self,$chrom,$line) = @_;

        chomp $line;
        my ($sub_line) = ($line =~ /^([^#]*)/);
        my ($seq,$source,$type,$start,$end) = split(/\t/,$sub_line);

	die "error: no seq, source, or type" if (not defined $seq or not defined $source or not defined $type);
	die "error: type=$type is not 'region'" if ($type ne 'region');

	my ($real_start, $real_end) = $self->_parse_coords($chrom,$start,$end,$line);

        return ($real_start,$real_end);
}

sub _parse_tsv_line
{
	my ($self,$line) = @_;

        chomp $line;
        my ($sub_line) = ($line =~ /^([^#]*)/);
        my ($chrom,$start,$end) = split(/\t/,$sub_line);

	die "error: chrom not valid" if (not defined $chrom);
	die "error: start not valid" if (not defined $start);
	die "error: end not valid" if (not defined $end);

	my ($real_start, $real_end) = $self->_parse_coords($chrom,$start,$end,$line);

        return ($chrom,$real_start,$real_end);
}

sub _remove_bad_pos
{
	my ($self,$core_chrom,$core_pos,$bad_pos_file) = @_;

        # open bad positions file
	my $bad_positions_parser = InvalidPositions->new;
	my ($bad_positions,$total_bad_positions) = $bad_positions_parser->read_invalid_positions($bad_pos_file);
	for my $chrom_pos (keys %$bad_positions)
        {
        	delete $core_pos->{$chrom_pos};
        }
}

sub _parse_core_positions
{
	my ($self,$core_file,$core_format) = @_;

	my %core_positions;

	my $core_contig_name = undef;
	if ($core_format eq 'gff')
	{
		$core_contig_name = basename($core_file, '.gff');
		die "error: $core_file does not have proper name" if (not defined $core_contig_name);
	}

        # open core positions file
        open(my $cfh, "<$core_file") or die "Could not open $core_file: $!";
        while(my $line = readline($cfh))
        {
		my ($start,$end);
		if ($core_format eq 'gff')
		{
                	($start,$end) = $self->_parse_gff_line($core_contig_name,$line);
		}
		elsif ($core_format eq 'tsv')
		{
                	($core_contig_name,$start,$end) = $self->_parse_tsv_line($line);
		}
		else
		{
			die "error: invalid core_format=$core_format";
		}

                next if (not defined $start or not defined $end);

                for (my $i = $start; $i <= $end; $i++)
                {
                        $core_positions{"${core_contig_name}_${i}"} = 1;
                }
        }
        close($cfh);

	return ($core_contig_name,\%core_positions);
}

1;
