
## Cython-ization

Unlike [NREL SPA](https://midcdmz.nrel.gov/spa/), [NREL SOLPOS](https://www.nrel.gov/grid/solar-resource/solpos.html) and [SPCTRAL2](https://rredc.nrel.gov/solar//models/spectral/) seem not to face the sort of licensing restrictions that caused [pvlib](https://pvlib-python.readthedocs.io/en/stable/) to not include the C source files in [their Cython](https://github.com/pvlib/pvlib-python/tree/master/pvlib/spa_c_files). 

Therefore I am including the (mostly) untouched files here, as downloaded from those links. 
* SOLPOS on 24-Feb-20
* SPCTRAL2 on X

### Necessary edits: 
* 'spectrl2_2.c': `#include` -> `spectral2_2.h` not `spectrl2.h` (L7)
  This was necessary in order to build without error.

Between the different versions (C, Fortran, Excel) of SPCTRAL2, they [name the program slightly differently](https://www.nrel.gov/grid/solar-resource/spectral.html). I stick with "spctral2" because, although lowercase, it is faithful to the official name.

### Other edits:

* Adding solar zenith angle (`zenref`) from SOLPOS as a SPCTRAL2 output, in order to check the value being used.


### Cython code:
* `cspctral2.pxd`: external declarations (~ the header file(s))
* `spctral2.pyx`: 

### Compiling
```sh
python setup.py build_ext --inplace
```
to compile for use in the current directory *only* (not added to PATH). This requires an available C compiler (easy to get on Linux, or on Mac with homebrew, [less clear on Windows](https://github.com/cython/cython/wiki/CythonExtensionsOnWindows)).

