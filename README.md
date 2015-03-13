vcf2pseudoalignment
===================

This project contains a number of PERL scripts for converting VCF files to an alignment of SNPs.  This is a subproject of the larger core phylogenomics pipeline project located at https://github.com/apetkau/core-phylogenomics.

Installation
============

Following the install guide of the parent project should automatically download this software.  However, if you wish to install this project independently please follow the below steps.

Step 1: Download Software
-------------------------

Download the software with:

	$ git clone https://github.com/apetkau/vcf2pseudoalignment.git

Step 2: Install Dependencies
----------------------------

This project requires [vcftools](http://vcftools.sourceforge.net/) as a dependency.  Please download and install this software.

In addition, this software requires the Perl modules from vcftools someplace within your library path (in PERL5LIB for example).  This can be accomplished by running, for example:

	$ export PERL5LIB=/path/to/vcftools/lib/perl5/site_perl/:$PERL5LIB

Or, you can simply link up the appropriate perl modules within the **vcf2pseudoalign/lib** directory.  For example:

	$ ln -s /path/to/vcftools/lib/perl5/site_perl/*.pm /path/to/vcf2pseudoalign/lib
	$ ls /path/to/vcf2pseudoalign/lib
	Align  CorePositions.pm  FaSlice.pm  InvalidPositions.pm  NucmerPositionsChecker.pm  PositionsTable.pm  Vcf.pm  VcfStats.pm

Need to compile custom bcftools plugin and modify a few ENV variables
	$ cd bcfplugins/bcftools-1.2
	$ make
	$ export PATH=`pwd`:$PATH
	$ export LD_LIBRARY_PATH=`pwd`/htslib-1.2.1:$LD_LIBRARY_PATH
	$ export BCFTOOLS_PLUGINS=`pwd`/plugins:$BCFTOOLS_PLUGINS

	
Step 3: Run Tests
-----------------

In order to run the tests, please run the command:

	$ cd /path/to/vcf2pseudoalign
	$ prove
	t/compare_pseudoalign_nucmer.t .. ok    
	t/extract_snps_metaalign.t ...... ok   
	t/find-positions-used.t ......... ok   
	t/nucmer_align.t ................ ok    
	t/variant_calls.t ............... ok     
	t/vcf2core.t .................... ok    
	All tests successful.
	Files=6, Tests=232, 41 wallclock secs ( 0.28 usr  0.02 sys + 30.71 cusr  7.60 csys = 38.61 CPU)
	Result: PASS

