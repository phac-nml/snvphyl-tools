SNVPhyl Tools
=============

This project contains a number of tools for the SNVPhyl whole genome phylogeny pipeline.

Installation
============

Please follow the below steps.

Step 1: Download Software
-------------------------

Download the software with:

```bash
git clone http://irida.corefacility.ca/gitlab/analysis-pipelines/snvphyl-tools.git
```

Step 2: Install Dependencies
----------------------------

* [vcftools](http://vcftools.sourceforge.net/)
* [bcftools](http://www.htslib.org/download/)
* [BioPerl](http://www.bioperl.org/wiki/Main_Page)

In addition, this software requires the Perl modules from vcftools someplace within your library path (in PERL5LIB for example).  This can be accomplished by running, for example:

	$ export PERL5LIB=/path/to/vcftools/lib/perl5/site_perl/:$PERL5LIB

Or, you can simply link up the appropriate perl modules within the **vcf2pseudoalign/lib** directory.  For example:

	$ ln -s /path/to/vcftools/lib/perl5/site_perl/*.pm /path/to/vcf2pseudoalign/lib
	$ ls /path/to/vcf2pseudoalign/lib
	Align  CorePositions.pm  FaSlice.pm  InvalidPositions.pm  NucmerPositionsChecker.pm  PositionsTable.pm  Vcf.pm  VcfStats.pm

Need to compile custom bcftools plugin and modify a few ENV variables
```
	cd bcfplugins/bcftools-1.3
	make
	export PATH=`pwd`:$PATH
	export LD_LIBRARY_PATH=`pwd`/htslib-1.3:$LD_LIBRARY_PATH
	export BCFTOOLS_PLUGINS=`pwd`/plugins:$BCFTOOLS_PLUGINS
```

Step 3: Run Tests
-----------------

In order to run the tests, please run the command:

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
