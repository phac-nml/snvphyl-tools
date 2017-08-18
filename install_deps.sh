#!/bin/bash

set -e


#flag to silent all make and configures being ran.
#set to empty string if you want to see everything in its glory
silent="--silent"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"



echo "Fetching and installing htslib...";

#installing all dependencies for building bcftools and our plugin
htslib='htslib-1.5'

if [ ! -d "$htslib" ]; then

    rm -rf $htslib
fi


curl -s https://codeload.github.com/samtools/htslib/tar.gz/1.5 > htslib.tar.gz
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

bcftools='bcftools-1.5'

if [ ! -d "$bcftools" ]; then

    rm -rf $bcftools
fi



curl -s https://codeload.github.com/samtools/bcftools/tar.gz/1.5 > bcftools.tar.gz
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

echo "Fetching and installing Mummer"

mummer='MUMmer3.23'

if [ ! -d "$mummer" ]; then

    rm -rf $mummer
fi


#fetching and putting on the path mummer for find-repeats.pl script
curl -s -L https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz/download > mummer.tar.gz
tar -zxf mummer.tar.gz
rm mummer.tar.gz

cd $mummer
#completely silenting the output from make
make $silent 2> /dev/null

#add mummer to the path
toolPATH=$toolPATH:`pwd`

cd $DIR

echo "Done!";
echo "Both PATH and BCFTOOL_PLUGINS have been modified for this terminal session."
echo "Please add following environment variables to your ~/.bashrc if long term use"
echo "export PATH=$toolPATH:\$PATH"
echo "export BCFTOOLS_PLUGINS=$plugins"


export PATH=$toolPATH:$PATH
export BCFTOOLS_PLUGINS=$plugins
