package NucmerPositionsChecker;

use warnings;
use strict;

use Set::Scalar;
use Align::Nucmer;

sub new
{
	my ($class,$reference_file,$genome_file,$genome_name,$pipeline_core_snp,$core_positions,$verbose) = @_;

	my $self = {};
	bless($self,$class);

	die "error: reference_file is not defined" if (not defined $reference_file);
	die "error: reference_file=$reference_file does not exist" if (not -e $reference_file);
	die "error: genome_file is not defined" if (not defined $genome_file);
	die "error: genome_file=$genome_file does not exist" if (not -e $genome_file);
	die "error: genome_name is not defined" if (not defined $genome_name);
	die "error: pipeline_core_snp is not defined" if (not defined $pipeline_core_snp);
	die "error: core_positions is not defined" if (not defined $core_positions);

	$verbose = 0 if (not defined $verbose);
	$self->{'verbose'} = $verbose;

	$self->_parse_genome_nucmer($reference_file, $genome_file, $genome_name, $pipeline_core_snp, $core_positions);

	return $self;
}

sub _parse_genome_nucmer
{
	my ($self, $reference_file, $genome_file, $genome_name, $pipeline_core_snp, $core_positions) = @_;
	my $verbose = $self->{'verbose'};

	my $nucmer_snp_set = Set::Scalar->new;
	my $nucmer_snp_set_core_pos = Set::Scalar->new;
	my $nucmer_ref_set = Set::Scalar->new;
	my $nucmer_ref_set_core_pos = Set::Scalar->new;

	my $nucmer_align_parser = Align::Nucmer->new($verbose);
	my $nucmer_results = $nucmer_align_parser->align_and_parse($reference_file,$genome_file);
	die "error: could not parse nucmer results for $reference_file vs. $genome_file" if (not defined $nucmer_results);

	my $nucmer_snps = $nucmer_results->{'snps'};
	$self->_insert_nucmer_snps($nucmer_snps,$nucmer_snp_set,$nucmer_snp_set_core_pos,$core_positions);

	my $pipeline_snps = $pipeline_core_snp->{$genome_name};
	$self->_insert_nucmer_reference($pipeline_snps,$nucmer_ref_set,$nucmer_ref_set_core_pos,$nucmer_snp_set,$nucmer_results,$core_positions);


	$self->{'nucmer_snp_set'} = $nucmer_snp_set;
	$self->{'nucmer_snp_set_core_pos'} = $nucmer_snp_set_core_pos;
	$self->{'nucmer_ref_set'} = $nucmer_ref_set;
	$self->{'nucmer_ref_set_core_pos'} = $nucmer_ref_set_core_pos;
}

sub get_nucmer_snp_set
{
	my ($self) = @_;

	return $self->{'nucmer_snp_set'};
}

sub get_nucmer_snp_set_core_pos
{
	my ($self) = @_;

	return $self->{'nucmer_snp_set_core_pos'};
}

sub get_nucmer_ref_set_core_pos
{
	my ($self) = @_;

	return $self->{'nucmer_ref_set_core_pos'};
}

sub get_nucmer_ref_set
{
	my ($self) = @_;

	return $self->{'nucmer_ref_set'};
}

sub _insert_nucmer_reference
{
	my ($self,$pipeline_snps,$nucmer_set,$nucmer_set_core_pos,$nucmer_snp_set,$nucmer_results,$core_positions) = @_;

	for my $chrom (keys %$pipeline_snps)
	{
		my $positions = $pipeline_snps->{$chrom};
		for my $pos (keys %$positions)
		{
			my $nucmer_results_all = $nucmer_results->{'all'};
			die "error: nucmer_results table is invalid" if (not defined $nucmer_results_all);
			my $nucmer_results_for_pos = $nucmer_results_all->{$chrom}->{$pos};
			my $nucmer_ref;
			my $nucmer_alt;
			my $pipeline_ref = $positions->{$pos}->{'reference'};
			my $pipeline_alt = $positions->{$pos}->{'alternative'};;
			die "error: no pipline_ref or pipeline_alt defined" if (not defined $pipeline_ref or not defined $pipeline_alt);

			if (not defined $nucmer_results_for_pos)
			{
				$nucmer_ref = 'UNKNOWN';
				$nucmer_alt = 'UNKNOWN';

				print STDERR "warning: no nucmer results for $chrom:$pos, using 'UNKNOWN'\n";
			}
			else
			{
				$nucmer_ref = $nucmer_results_for_pos->{'ref'};
				$nucmer_alt = $nucmer_results_for_pos->{'alt'};
			}

			# if this was a snp, should have been found in the nucmer snp results
			if ($nucmer_ref ne $nucmer_alt)
			{
				if (not $nucmer_snp_set->has("$chrom\t$pos\t$nucmer_ref\t$nucmer_alt"))
				{
					die "error: snp $chrom:$pos:$nucmer_ref:$nucmer_alt found from nucmer align results not present in show-snps results";
				}
			}
			else
			{
				$nucmer_set->insert("$chrom\t$pos\t$nucmer_ref\t$nucmer_alt");
				if (exists $core_positions->{"${chrom}_$pos"})
				{
					$nucmer_set_core_pos->insert("$chrom\t$pos\t$nucmer_ref\t$nucmer_alt");
				}
			}
		}
	}
}

sub _insert_nucmer_snps
{
	my ($self,$nucmer_snps,$nucmer_set,$nucmer_set_core_pos,$core_positions) = @_;

	my $verbose = $self->{'verbose'};

	for my $contig (keys %$nucmer_snps)
	{
		my $positions = $nucmer_snps->{$contig};
		for my $pos (keys %$positions)
		{
			my $snp_data  = $positions->{$pos};
			my $ref = $snp_data->{'ref'};
			my $alt = $snp_data->{'alt'};

			if ($ref eq '.' or $alt eq '.')
			{
				print STDERR "Indel found at $contig:$pos:ref=$ref:alt=$alt, skipping\n" if ($verbose);
			}
			else
			{
				$nucmer_set->insert("$contig\t$pos\t$ref\t$alt");
				if (exists $core_positions->{"${contig}_$pos"})
				{
					$nucmer_set_core_pos->insert("$contig\t$pos\t$ref\t$alt");
				}
			}
		}
	}
}

1;
