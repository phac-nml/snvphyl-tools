package NucmerPositionsChecker;

use warnings;
use strict;

use Set::Scalar;
use Align::Nucmer;

sub new
{
	my ($class,$reference_file,$genome_file,$genome_name,$pipeline_core_snv,$core_positions,$verbose) = @_;

	my $self = {};
	bless($self,$class);

	die "error: reference_file is not defined" if (not defined $reference_file);
	die "error: reference_file=$reference_file does not exist" if (not -e $reference_file);
	die "error: genome_file is not defined" if (not defined $genome_file);
	die "error: genome_file=$genome_file does not exist" if (not -e $genome_file);
	die "error: genome_name is not defined" if (not defined $genome_name);
	die "error: pipeline_core_snv is not defined" if (not defined $pipeline_core_snv);
	die "error: core_positions is not defined" if (not defined $core_positions);

	$verbose = 0 if (not defined $verbose);
	$self->{'verbose'} = $verbose;

	$self->_parse_genome_nucmer($reference_file, $genome_file, $genome_name, $pipeline_core_snv, $core_positions);

	return $self;
}

sub _parse_genome_nucmer
{
	my ($self, $reference_file, $genome_file, $genome_name, $pipeline_core_snv, $core_positions) = @_;
	my $verbose = $self->{'verbose'};

	my $nucmer_snv_set = Set::Scalar->new;
	my $nucmer_snv_set_core_pos = Set::Scalar->new;
	my $nucmer_ref_set = Set::Scalar->new;
	my $nucmer_ref_set_core_pos = Set::Scalar->new;

	my $nucmer_align_parser = Align::Nucmer->new($verbose);
	my $nucmer_results = $nucmer_align_parser->align_and_parse($reference_file,$genome_file);
	die "error: could not parse nucmer results for $reference_file vs. $genome_file" if (not defined $nucmer_results);

	my $nucmer_snvs = $nucmer_results->{'snvs'};
	$self->_insert_nucmer_snvs($nucmer_snvs,$nucmer_snv_set,$nucmer_snv_set_core_pos,$core_positions);

	my $pipeline_snvs = $pipeline_core_snv->{$genome_name};
	$self->_insert_nucmer_reference($pipeline_snvs,$nucmer_ref_set,$nucmer_ref_set_core_pos,$nucmer_snv_set,$nucmer_snv_set_core_pos,$nucmer_results,$core_positions);


	$self->{'nucmer_snv_set'} = $nucmer_snv_set;
	$self->{'nucmer_snv_set_core_pos'} = $nucmer_snv_set_core_pos;
	$self->{'nucmer_ref_set'} = $nucmer_ref_set;
	$self->{'nucmer_ref_set_core_pos'} = $nucmer_ref_set_core_pos;
}

sub get_nucmer_snv_set
{
	my ($self) = @_;

	return $self->{'nucmer_snv_set'};
}

sub get_nucmer_snv_set_core_pos
{
	my ($self) = @_;

	return $self->{'nucmer_snv_set_core_pos'};
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

sub get_unknown_pipeline_positions
{
	my ($self) = @_;

	return $self->{'unknown_pipeline_positions'};
}

sub get_unknown_snvs_removed_count
{
	my ($self) = @_;

	return $self->{'unknown_snvs_removed_count'};
}

sub _insert_nucmer_reference
{
	my ($self,$pipeline_snvs,$nucmer_set,$nucmer_set_core_pos,$nucmer_snv_set,$nucmer_snv_set_core_pos,$nucmer_results,$core_positions) = @_;

	my $unknown_pipeline_positions = Set::Scalar->new;
	my $unknown_snvs_removed_count = 0; # snvs removed when removing unknown positions

	for my $chrom (keys %$pipeline_snvs)
	{
		my $positions = $pipeline_snvs->{$chrom};
		for my $pos (keys %$positions)
		{
			my $nucmer_results_all = $nucmer_results->{'all'};
			die "error: nucmer_results table is invalid" if (not defined $nucmer_results_all);
			my $nucmer_results_for_pos = $nucmer_results_all->{$chrom}->{$pos};
			my $pipeline_ref = $positions->{$pos}->{'reference'};
			my $pipeline_alt = $positions->{$pos}->{'alternative'};;
			die "error: no pipline_ref or pipeline_alt defined" if (not defined $pipeline_ref or not defined $pipeline_alt);

			if (not defined $nucmer_results_for_pos)
			{
				print STDERR "warning: no nucmer results for $chrom:$pos, adding to uknown set\n";
				$unknown_pipeline_positions->insert("$chrom\t$pos\t$pipeline_ref\t$pipeline_alt");
				$unknown_snvs_removed_count++ if ($pipeline_ref ne $pipeline_alt);
			}
			else
			{
				my $nucmer_ref = $nucmer_results_for_pos->{'ref'};
				my $nucmer_alt = $nucmer_results_for_pos->{'alt'};

				# if this was a snv, should have been found in the nucmer snv results
				if ($nucmer_ref ne $nucmer_alt)
				{
					if (not $nucmer_snv_set->has("$chrom\t$pos\t$nucmer_ref\t$nucmer_alt"))
					{
						print STDERR "warning: snv $chrom:$pos:$nucmer_ref:$nucmer_alt found from nucmer align results not present in show-snvs results\n";
						$nucmer_snv_set->insert("$chrom\t$pos\t$nucmer_ref\t$nucmer_alt");
						if (exists $core_positions->{"${chrom}_$pos"})
						{
							$nucmer_snv_set_core_pos->insert("$chrom\t$pos\t$nucmer_ref\t$nucmer_alt");
						}
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

	$self->{'unknown_pipeline_positions'} = $unknown_pipeline_positions;
	$self->{'unknown_snvs_removed_count'} = $unknown_snvs_removed_count;
}

sub _insert_nucmer_snvs
{
	my ($self,$nucmer_snvs,$nucmer_set,$nucmer_set_core_pos,$core_positions) = @_;

	my $verbose = $self->{'verbose'};

	for my $contig (keys %$nucmer_snvs)
	{
		my $positions = $nucmer_snvs->{$contig};
		for my $pos (keys %$positions)
		{
			my $snv_data  = $positions->{$pos};
			my $ref = $snv_data->{'ref'};
			my $alt = $snv_data->{'alt'};

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
