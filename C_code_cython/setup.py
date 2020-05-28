from distutils.core import setup
from distutils.extension import Extension
import os

from Cython.Build import cythonize

DIRNAME = os.path.dirname(__file__)

# sources = [os.path.join(DIRNAME, src) for src in ('spctral2.pyx', 'spectrl2_2.c')]
sources = [os.path.join(DIRNAME, src) for src in ('spctral2.pyx', 'spectrl2_2.c', 'solpos.c')]

# create specify extension module name using Extension
extensions = [Extension('csp2_py', sources)]

# cythonize
extensions = cythonize(extensions, language_level=3)

setup(
    name="Cython interface to original C version of NREL's SPCTRAL2", 
    ext_modules=extensions,
)
