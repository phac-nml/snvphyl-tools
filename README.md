SNVPhyl Tools
=============

This project will download and compile a number of dependency tools for the [SNVPhyl][] whole genome phylogeny pipeline.

Legal
=====

Copyright 2012-2017 Government of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.


## Install Options


### 1. Install with bioconda (recommended)

```shell
conda install snvphyl-tools
```


### 2. Install dependencies using install_deps.sh
Script will download and compile in the current working directory [bcftools](http://www.htslib.org/),[htslib](http://www.htslib.org/), [Mummer](http://mummer.sourceforge.net/) and [Vcftools](https://vcftools.github.io/index.html).


#### Step a: Install Dependencies
----------------------------



```bash
git clone https://github.com/phac-nml/snvphyl-tools.git
cd snvphyl-tools
./install_deps.sh
#wait patiently to completion and setup the following environment variables
#example below
export PATH=/home/test/snvphyl-tools/bcftools-1.9:/home/test/snvphyl-tools/MUMmer3.23:$PATH
export BCFTOOLS_PLUGINS=/home/test//snvphyl-tools/bcftools-1.9/plugins
PERL5LIB=/home/test/snvphyl-tools/lib/perl5
```

#### Step b: Run Tests
-----------------

In order to run the tests, please run the command:

	$ prove
	t/extract_snvs_metaalign.t ...... ok   
	t/find-positions-used.t ......... ok   
	t/variant_calls.t ............... ok     
	t/vcf2core.t .................... ok    
	All tests successful.
	Files=6, Tests=232, 41 wallclock secs ( 0.28 usr  0.02 sys + 30.71 cusr  7.60 csys = 38.61 CPU)
	Result: PASS


[SNVPhyl]: http://snvphyl.readthedocs.org/
[filter_snv_density]: bcfplugins/filter_snv_density.c
