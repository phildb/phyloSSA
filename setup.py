#!/Users/philippe/anaconda/bin/python
from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "PhyloSSA2",
    ext_modules = cythonize('phyloSSA2.py'),  # accepts a glob pattern
)
