"""
Run the SPCTRAL2 C program directly with Python 
with the help of Cython.

This module defines the runner function.
"""
# cimport cython
# cimport numpy as cnp
# import numpy as np

cimport cspctral2


def run_spctral2(*, 
    lat, lon,
    year, month, day, hour, minute, second,
    utcoffset,
    temp, press, 
    tau500, watvap, ozone
    ):
    """Run SPCTRAL2 C version.

    This is essentially a copy of the 'spectest.c' program, 
    except that we can keep the return values in memory instead of writing them out.

    When run in Python, this function returns a `dict`, where the keys correspond to 
    the `struct` members defined in `cspctral2.pxd`
    """

    cdef cspctral2.specdattype sp2
    # ctypedef struct cspctral2.specdattype sp2
 
    # cdef cspctral2.specdattype sp2, *specdat

    # specdat = &sp2; 

    cspctral2.S_spec_init(&sp2)
    # cspctral2.S_spec_init(specdat)

    # option 1: spectral irradiance (W/m^2/um) and wavelength in (um)
    sp2.units = 1

    sp2.latitude = lat
    sp2.longitude = lon

    sp2.year = year
    sp2.month = month
    sp2.day = day
    sp2.hour = hour
    sp2.minute = minute
    sp2.second = second

    sp2.timezone = utcoffset

    sp2.temp = temp
    sp2.press = press

    sp2.tau500 = tau500
    sp2.watvap = watvap
    sp2.ozone = ozone

    # this should be parallel to (spherical Earth) ground for NH at least
    sp2.tilt = lat
    sp2.aspect = 0.  # (S=180, N=0, E=90, W=270); aspect (azimuth angle) of collector panel

    err = cspctral2.S_spectral2(&sp2)
    # err = cspctral2.S_spectral2(specdat)

    return sp2
