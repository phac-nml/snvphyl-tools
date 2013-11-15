#!/usr/bin/env perl

use warnings;
use strict;

use Test::More;
use Test::Deep;
use File::Temp 'tempfile';

use FindBin;
use lib $FindBin::Bin.'/../lib';

use Align::Nucmer;

sub build_in_temp
{
	my ($string) = @_;

	my ($fh,$file) = tempfile('align_nucmer.XXXXXX',TMPDIR => 1, UNLINK => 1);
	print $fh $string;
	close($fh);

	return $file;
}

sub test_parse_alignments_1
{
	my $show_aligns =
'reference.fasta query.fasta

============================================================
-- Alignments between ref and query

-- BEGIN alignment [ +1 1 - 95 | +1 1 - 95 ]


1          tgaaatcgaatcggattcgaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
1          tgaattcgaatcggattcgaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
               ^                                            

50         aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
50         aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
                                                         


--   END alignment [ +1 1 - 95 | +1 1 - 95 ]

============================================================';

	my $show_snps =
'reference.fasta query.fasta
NUCMER

[P1]	[SUB]	[SUB]	[P2]	[BUFF]	[DIST]	[FRM]	[TAGS]
5	A	T	5	5	5	1	1	ref	query';
	my $show_aligns_file = build_in_temp($show_aligns);
	my $show_snps_file = build_in_temp($show_snps);

	my $expected_out = {
		'all' =>
			{ 'ref' =>
				{
					1 => { "ref" => "T", "alt" => "T" },
					2 => { "ref" => "G", "alt" => "G" },
					3 => { "ref" => "A", "alt" => "A" },
					4 => { "ref" => "A", "alt" => "A" },
					5 => { "ref" => "A", "alt" => "T" },
					6 => { "ref" => "T", "alt" => "T" },
					7 => { "ref" => "C", "alt" => "C" },
					8 => { "ref" => "G", "alt" => "G" },
					9 => { "ref" => "A", "alt" => "A" },
					10 => { "ref" => "A", "alt" => "A" },
					11 => { "ref" => "T", "alt" => "T" },
					12 => { "ref" => "C", "alt" => "C" },
					13 => { "ref" => "G", "alt" => "G" },
					14 => { "ref" => "G", "alt" => "G" },
					15 => { "ref" => "A", "alt" => "A" },
					16 => { "ref" => "T", "alt" => "T" },
					17 => { "ref" => "T", "alt" => "T" },
					18 => { "ref" => "C", "alt" => "C" },
					19 => { "ref" => "G", "alt" => "G" },
					20 => { "ref" => "A", "alt" => "A" },
					21 => { "ref" => "A", "alt" => "A" },
					22 => { "ref" => "A", "alt" => "A" },
					23 => { "ref" => "A", "alt" => "A" },
					24 => { "ref" => "A", "alt" => "A" },
					25 => { "ref" => "A", "alt" => "A" },
					26 => { "ref" => "A", "alt" => "A" },
					27 => { "ref" => "A", "alt" => "A" },
					28 => { "ref" => "A", "alt" => "A" },
					29 => { "ref" => "A", "alt" => "A" },
					30 => { "ref" => "A", "alt" => "A" },
					31 => { "ref" => "A", "alt" => "A" },
					32 => { "ref" => "A", "alt" => "A" },
					33 => { "ref" => "A", "alt" => "A" },
					34 => { "ref" => "A", "alt" => "A" },
					35 => { "ref" => "A", "alt" => "A" },
					36 => { "ref" => "A", "alt" => "A" },
					37 => { "ref" => "A", "alt" => "A" },
					38 => { "ref" => "A", "alt" => "A" },
					39 => { "ref" => "A", "alt" => "A" },
					40 => { "ref" => "A", "alt" => "A" },
					41 => { "ref" => "A", "alt" => "A" },
					42 => { "ref" => "A", "alt" => "A" },
					43 => { "ref" => "A", "alt" => "A" },
					44 => { "ref" => "A", "alt" => "A" },
					45 => { "ref" => "A", "alt" => "A" },
					46 => { "ref" => "A", "alt" => "A" },
					47 => { "ref" => "A", "alt" => "A" },
					48 => { "ref" => "A", "alt" => "A" },
					49 => { "ref" => "A", "alt" => "A" },
					50 => { "ref" => "A", "alt" => "A" },
					51 => { "ref" => "A", "alt" => "A" },
					52 => { "ref" => "A", "alt" => "A" },
					53 => { "ref" => "A", "alt" => "A" },
					54 => { "ref" => "A", "alt" => "A" },
					55 => { "ref" => "A", "alt" => "A" },
					56 => { "ref" => "A", "alt" => "A" },
					57 => { "ref" => "A", "alt" => "A" },
					58 => { "ref" => "A", "alt" => "A" },
					59 => { "ref" => "A", "alt" => "A" },
					60 => { "ref" => "A", "alt" => "A" },
					61 => { "ref" => "A", "alt" => "A" },
					62 => { "ref" => "A", "alt" => "A" },
					63 => { "ref" => "A", "alt" => "A" },
					64 => { "ref" => "A", "alt" => "A" },
					65 => { "ref" => "A", "alt" => "A" },
					66 => { "ref" => "A", "alt" => "A" },
					67 => { "ref" => "A", "alt" => "A" },
					68 => { "ref" => "A", "alt" => "A" },
					69 => { "ref" => "A", "alt" => "A" },
					70 => { "ref" => "A", "alt" => "A" },
					71 => { "ref" => "A", "alt" => "A" },
					72 => { "ref" => "A", "alt" => "A" },
					73 => { "ref" => "A", "alt" => "A" },
					74 => { "ref" => "A", "alt" => "A" },
					75 => { "ref" => "A", "alt" => "A" },
					76 => { "ref" => "A", "alt" => "A" },
					77 => { "ref" => "A", "alt" => "A" },
					78 => { "ref" => "A", "alt" => "A" },
					79 => { "ref" => "A", "alt" => "A" },
					80 => { "ref" => "A", "alt" => "A" },
					81 => { "ref" => "A", "alt" => "A" },
					82 => { "ref" => "A", "alt" => "A" },
					83 => { "ref" => "A", "alt" => "A" },
					84 => { "ref" => "A", "alt" => "A" },
					85 => { "ref" => "A", "alt" => "A" },
					86 => { "ref" => "A", "alt" => "A" },
					87 => { "ref" => "A", "alt" => "A" },
					88 => { "ref" => "A", "alt" => "A" },
					89 => { "ref" => "A", "alt" => "A" },
					90 => { "ref" => "A", "alt" => "A" },
					91 => { "ref" => "A", "alt" => "A" },
					92 => { "ref" => "A", "alt" => "A" },
					93 => { "ref" => "A", "alt" => "A" },
					94 => { "ref" => "A", "alt" => "A" },
					95 => { "ref" => "A", "alt" => "A" }
				}
			},
		'snps' =>
			{ 'ref' =>
				{ 5 => {'ref' => 'A', 'alt' => 'T'}}
			}
		};
		
	my $align_nucmer = Align::Nucmer->new;
	my $actual_out = $align_nucmer->parse_alignments($show_aligns_file,$show_snps_file);

	cmp_deeply($actual_out,$expected_out,'test_parse_alignments_1 output passes');
}

### MAIN ###
print "Testing Align::Nucmer package\n";
test_parse_alignments_1;

done_testing();
