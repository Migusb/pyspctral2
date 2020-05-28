"""
External declarations from C header file `spectrl2_2.h`

"Redefinition of the C API in a .pxd file"
  https://cython.readthedocs.io/en/latest/src/tutorial/clibraries.html#defining-external-declarations

Notes
-----
To run spctral2 (as seen in file 'spectest.c'), requires only the headers in 'spectrl2.h'
however, 'spectrl2_2.c' requires both 'spectrl2.h' and 'solpos00.h'
"""

cdef extern from "spectrl2_2.h":

    #> the struct that collects all of the input parameters
    # ctypedef struct specdattype:
    cdef struct specdattype:
        
        int units

        int year
        int month
        int day
        int hour
        int minute
        int second
        float timezone

        float alpha
        
        float tilt
        float aspect

        float assym

        float latitude
        float longitude

        float temp
        float press

        float[122] specdif
        float[122] specdir
        float[122] specetr
        float[122] specglo
        float[122] specx

        float[6] spcwvr
        
        float tau500
        float watvap
        float ozone

        
    # void* specdat

    #> func to run the model for the conditions outlined in the specdat struct
    extern int S_spectral2(specdattype *specdat)
    
    #> func to initialize struct
    void S_spec_init(specdattype *specdat)
