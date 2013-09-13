package CorePositions;

use warnings;
use strict;

use File::Basename;

sub new
{
	my ($class, $core_gff_file, $bad_pos_file) = @_;

	my $self = {};
	bless($self,$class);

	die "error: core_gff_file is not defined" if (not defined $core_gff_file);
	die "error: $core_gff_file does not exist" if (not -e $core_gff_file);
	die "error: $core_gff_file not a proper gff file" if ($core_gff_file !~ /\.gff$/i);

	die "error: bad_pos_file is not defined" if (not defined $bad_pos_file);
	die "error: $bad_pos_file does not exist" if (not -e $bad_pos_file);

	my ($core_chrom,$core_pos) = $self->_parse_core_positions($core_gff_file);
	die "error: core chromosome is undefined" if (not defined $core_chrom);

	$self->{'core_pos_count'} = scalar(keys %$core_pos);

	$self->_remove_bad_pos($core_chrom,$core_pos,$bad_pos_file);

	$self->{'used_pos_count'} = scalar(keys %$core_pos);

	$self->{'core_gff_file'} = $core_gff_file;
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

sub _parse_bad_pos_line
{
	my ($self,$line) = @_;

        chomp $line;
        my ($sub_line) = ($line =~ /^([^#]*)/);
        my ($chrom,$start,$end) = split(/\t/,$sub_line);

	my ($real_start, $real_end) = $self->_parse_coords($chrom,$start,$end,$line);

        return ($chrom,$real_start,$real_end);
}

sub _remove_bad_pos
{
	my ($self,$core_chrom,$core_pos,$bad_pos_file) = @_;

        # open bad positions file
        open(my $bfh, "<$bad_pos_file") or die "Could not open $bad_pos_file: $!";
        while(my $line = readline($bfh))
        {
                my ($chrom,$start,$end) = $self->_parse_bad_pos_line($line);
                next if (not defined $chrom or not defined $start or not defined $end);

		die "error: chromosome name: $core_chrom from GFF file not same as".
		" chromosome name: $chrom from bad_pos file" if ($core_chrom ne $chrom);

                # remove any bad positions from core
                for (my $i = $start; $i <= $end; $i++)
                {
                        delete $core_pos->{"${chrom}_${i}"}
                                if (exists $core_pos->{"${chrom}_${i}"});
                }
        }
        close($bfh);
}

sub _parse_core_positions
{
	my ($self,$core_gff_file) = @_;

	my %core_positions;

	my $core_contig_name = basename($core_gff_file, '.gff');
	die "error: $core_gff_file does not have proper name" if (not defined $core_contig_name);

        # open core positions file
        open(my $cfh, "<$core_gff_file") or die "Could not open $core_gff_file: $!";
        while(my $line = readline($cfh))
        {
                my ($start,$end) = $self->_parse_gff_line($core_contig_name,$line);
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
