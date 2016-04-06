SNVPhyl Tools
=============

This project contains a number of dependency tools for the [SNVPhyl][] whole genome phylogeny pipeline.

Legal
=====

Copyright 2012-2016 Government of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

Installation
============

Please follow the below steps.

Step 1: Download Software
-------------------------

Download the software with:

```bash
git clone https://irida.corefacility.ca/analysis-pipelines/snvphyl-tools.git
```

Step 2: Install Dependencies
----------------------------

* [vcftools](http://vcftools.sourceforge.net/)
* [bcftools](http://www.htslib.org/download/)
* [BioPerl](http://www.bioperl.org/wiki/Main_Page)

In addition, this software requires the Perl modules from vcftools someplace within your library path (in PERL5LIB for example).  This can be accomplished by running, for example:

	$ export PERL5LIB=/path/to/vcftools/lib/perl5/site_perl/:$PERL5LIB

Or, you can simply link up the appropriate perl modules within the **snvphyl-tools/lib/** directory.  For example:

	$ ln -s /path/to/vcftools/lib/perl5/site_perl/*.pm /path/to/snvphyl-tools/lib
	$ ls /path/to/snvphyl-tools/lib
	Align  CorePositions.pm  FaSlice.pm  InvalidPositions.pm  NucmerPositionsChecker.pm  PositionsTable.pm  Vcf.pm  VcfStats.pm

You will need to compile custom bcftools plugin and modify a few ENV variables.
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
	t/compare_snv_align_nucmer.t .. ok    
	t/extract_snvs_metaalign.t ...... ok   
	t/find-positions-used.t ......... ok   
	t/nucmer_align.t ................ ok    
	t/variant_calls.t ............... ok     
	t/vcf2core.t .................... ok    
	All tests successful.
	Files=6, Tests=232, 41 wallclock secs ( 0.28 usr  0.02 sys + 30.71 cusr  7.60 csys = 38.61 CPU)
	Result: PASS


[SNVPhyl]: http://snvphyl.readthedocs.org/
