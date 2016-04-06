package PositionsTable;

use warnings;
use strict;

# builds a new PositionsTable parser object
# input
#	verbose:  Print verbose information.
sub new
{
	my ($class,$verbose) = @_;

	my $self = {};
	$verbose = 0 if (not defined $verbose);
	$self->{'_verbose'} = $verbose;
	bless($self,$class);

	return $self;
}

# parses the PositionsTable and returns a hash table of the data
# input
#	positions_file:  The positions table file
# output
# 	a genomes_core_snv table mapping
# 		strain_id => {chrom => {pos => {'reference' => $ref_base, 'alternative' => $alt_base}}}
#	a genomes_core_snv_count table mapping
#		strain_id => count_snvs_from_reference
sub read_table
{
        my ($self,$positions_file) = @_;
	my $verbose = $self->{'_verbose'};

        open(my $fh, "<$positions_file") or die "Could not open $positions_file: $!";

        my $line = readline($fh);
        chomp($line);

        die "Error: no header line defined in $positions_file" if ($line !~ /^#Chromosome\tPosition\tStatus\tReference/);
        my (undef,undef,undef,@strains) = split(/\t/,$line);
        die "Error: no strains defined in $positions_file" if (@strains <= 0);
        die "Error: reference not in correct column" if ($strains[0] ne 'Reference');
        my %genomes_core_snv;
        my %genomes_core_snv_count;

        # initialize empty table for each strain
        for my $strain (@strains)
        {
                $genomes_core_snv{$strain} = undef;
                $genomes_core_snv_count{$strain} = 0;
        }

	my $valid = 0;
        my $total = 0;
        while(my $line = readline($fh))
        {
                chomp $line;
                my @values = split(/\t/,$line);

                my ($chrom,$pos,$status,@dna) = @values;

                if (scalar(@dna) != scalar(@strains))
                {
                        die "Error: line $line does not have same number of entries as header for $positions_file";
                }
                elsif ($status ne 'valid')
                {
                        print STDERR "skipping over line \"$line\": invalid\n" if ($verbose);
                }
                else
                {
                        for (my $i = 1; $i < @dna; $i++)
                        {
                                $genomes_core_snv{$strains[$i]}{$chrom}{$pos} = {'reference' => $dna[0],
                                        'alternative' => $dna[$i]};

                                if ($dna[0] ne $dna[$i])
                                {
                                        $genomes_core_snv_count{$strains[$i]}++;
                                }

                        }
                        $valid++;
                }
                $total++;
        }
        close $fh;

        print STDERR "Kept $valid valid positions out of $total total positions\n" if ($verbose);

        return (\%genomes_core_snv,\%genomes_core_snv_count);
}

1;
