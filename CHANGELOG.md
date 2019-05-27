# Changes

## 1.8.2

*  Update bcftools from 1.5 to version 1.9
*  Update consolidate_vcf.pl to new standard of bcftools version 1.9

## 1.8.1

* Re-structured repository to use a `Makefile.PL` for handling installation/dependency checking. For integration with Bioconda.
* Added `install_deps.sh` script for installing dependencies.
* Removed copy of `bcftools` code in this repository for compiling plugins and instead use official distribution.
* Added `--out_strains` option to `verify_mapping_quality.pl` to produce a list of only the strain/sample names which did not pass the minimum coverage checks.
