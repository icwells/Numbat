# Numbat is a package of tools for k-mer based local alignment of metagenomic DNA-Seq data and for identifying viral hosts via codon frequency analysis.
## Please note that I no longer work in the lab I made this for, so this package should be considered unsupported.
### However, some of the included scripts may prove useful, so feel free to use them.

Copyright 2017 by Shawn Rupp

## Installation
Download the repository:

git clone https://github.com/icwells/Numbat.git

Most of the scripts are written in python3, but several contain Cython modules which
must be compiled. Cython can be installed from the pypi repository or via Miniconda 
(it is installed by default with the full Anaconda package).

### To install with Miniconda:
conda install cython

### Compiling Numbat:
cd Numbat/

./install.sh

## Please refer to NumbatReadMe.pdf for more detailed instructions on running the program
