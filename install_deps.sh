#!/bin/bash

set -e


#flag to silent all make and configures being ran.
#set to empty string if you want to see everything in its glory
silent="--silent"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"



echo "Fetching and installing htslib...";

#installing all dependencies for building bcftools and our plugin
htslib='htslib-1.9'

if [ ! -d "$htslib" ]; then

    rm -rf $htslib
fi


curl -s https://codeload.github.com/samtools/htslib/tar.gz/1.9 > htslib.tar.gz
tar -zxf htslib.tar.gz
#clean up
rm htslib.tar.gz

cd $htslib
autoreconf 
./configure --disable-bz2 --disable-lzma $silent
make $silent 2> /dev/null

export HTSLIB=`pwd`/htslib
cd $DIR

echo "Done!"

echo "Fetching and install bcftools"

bcftools='bcftools-1.9'

if [ ! -d "$bcftools" ]; then

    rm -rf $bcftools
fi



curl -s https://codeload.github.com/samtools/bcftools/tar.gz/1.9 > bcftools.tar.gz
tar -zxf bcftools.tar.gz

#clean up
rm bcftools.tar.gz

cp bcfplugins/filter_snv_density.c $bcftools/plugins
cd $bcftools

autoreconf
./configure $silent 2> /dev/null
make $silent 2> /dev/null

#add bcftools to the path
toolPATH=`pwd`

plugins=`pwd`/plugins

cd $DIR

echo "Done!"


echo "Installing cpan and other dependencies"

conda --version 1>/dev/null 2>/dev/null
exit_code = $?
if [ exit_code -ne 0 ]
	    echo "Please install 'conda' before continuing."
fi
#creating and activating test environment
conda create -y -n snvphyltesting
conda activate snvphyltesting

#creating new conda env with all perl dependencies and mummer,samtools and vcftools
conda install -y perl perl-bioperl perl-hash-merge perl-list-moreutils perl-math-round perl-parallel-forkmanager perl-string-util perl-template-toolkit perl-test-exception perl-text-csv perl-text-diff perl-vcftools-vcf mummer samtools vcftools perl-json-parse perl-string-util perl-parallel-forkmanager


echo "Please add following environment variables to the terminal or your ~/.bashrc for long term use"
echo "export PATH=$toolPATH:\$PATH"
echo "export BCFTOOLS_PLUGINS=$plugins"
export PATH=$toolPATH:$PATH
export BCFTOOLS_PLUGINS=$plugins
