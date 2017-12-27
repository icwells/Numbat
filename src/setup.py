'''This script will cythonize the Numbat package.'''

from distutils.core import setup
from Cython.Build import cythonize

FIO = "fastaIO.pyx"
KC = "kclust.pyx"

# Print blank lines to split output
print(("\n\tComipiling {}...\n").format(FIO))
setup(ext_modules=cythonize(FIO))
print(("\n\tComipiling {}...\n").format(KC))
setup(ext_modules=cythonize(KC))
print()
