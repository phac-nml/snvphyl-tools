use 5.006;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'SNVPhyl',
    EXE_FILES           => [ 'bin/snv_matrix.pl' ,
                             'bin/consolidate_vcfs.pl',
                             'bin/core-only.pl',
                             'bin/extract_snvs_metaalign.pl',
                             'bin/filter-positions.pl',
                             'bin/filter-stats.pl',
                             'bin/filter_unique_basepairs.pl',
                             'bin/filterVcf.pl',
                             'bin/find-positions-used.pl',
                             'bin/find-repeats.pl',
                             'bin/positions2phyloviz.pl',
                             'bin/positions2snv_alignment.pl',
                             'bin/positions2snv_invariant_alignment.pl',
                             'bin/rearrange_snv_matrix.pl',
                             'bin/ref_stats.pl',
                             'bin/reporter.pl',
                             'bin/snv_matrix.pl',
                             'bin/tsvToVcf.pl',
                             'bin/vcf2snv_alignment.pl',
                             'bin/verify_excess_coverage.pl',
                             'bin/verify_low_depth.pl',
                             'bin/verify_mapping_quality.pl'
                         ],    
    LICENSE          => 'apache',
    PL_FILES         => {},
    MIN_PERL_VERSION => '5.006',
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => '0',
    },
    BUILD_REQUIRES => {
        'Test::More' => '0',
    },
    PREREQ_PM => {
        'Bio::Perl' => '1.007001',
        'Switch' => '2.17',
        'List::MoreUtils' => '0.419',
        'JSON' => '2.94',
        'Test::JSON' => '0.11',
        'Math::Round' => '0.07',
        'Text::CSV' => '1.95',
        'Parallel::ForkManager' => '1.19',
        'Readonly' => '2.05',
        'Hash::Merge' => '0.200',
        'String::Util' => '1.26',
        'Bio::Phylo' => '0.58',
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'snvphyl-*' },
);
