#!/bin/bash

set -e


#flag to silent all make and configures being ran.
#set to empty string if you want to see everything in its glory
silent="--silent"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


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

#creating new conda env with all perl dependencies and mummer,samtools and vcftools
conda install -y perl perl-bioperl perl-hash-merge perl-list-moreutils perl-math-round perl-parallel-forkmanager perl-string-util perl-template-toolkit perl-test-exception perl-text-csv perl-text-diff perl-vcftools-vcf mummer samtools vcftools perl-json-parse perl-string-util perl-parallel-forkmanager


echo "Please add following environment variables to the terminal or your ~/.bashrc for long term use"
echo "export PATH=$toolPATH:\$PATH"
echo "export BCFTOOLS_PLUGINS=$plugins"
export PATH=$toolPATH:$PATH
export BCFTOOLS_PLUGINS=$plugins
