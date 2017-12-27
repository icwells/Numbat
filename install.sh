#!/bin/bash

##############################################################################
# This script will cythonize scripts for the Numbat package.
# 
# Required programs:	Cython
##############################################################################

FIO="fastaIO"
KC="kclust"

# Compile cython scripts and remove build files
cd src/
python setup.py build_ext --inplace
cd ../

mv src/$FIO.*.so bin/$FIO.so
mv src/$KC.*.so bin/$KC.so

rm -r src/build/
rm src/*.c
