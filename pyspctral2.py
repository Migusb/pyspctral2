
"""

Python module to run the SPCTRAL2 C code and catch output

"""

import datetime as dt
from glob import glob
import os
from subprocess import call

import numpy as np


C_template=r"""
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "spectrl2.h"   /* <-- include for spectrl2.  */


int main ( )
{{
        struct specdattype spdat, *specdat;  /* spectral2 data structure */
        
        /* point to the specral2 structure */
        specdat = &spdat;
        
        /* Initialize the data structure -- you can omit if ALL structure vaiables 
                        are assigned elswhere. Defaults in initialization can be overridden by 
                        reassigning below. */
           
        S_spec_init (specdat);
    
    /*  */
    
    specdat->longitude = {lon:.4f};  /* Note that latitude and longitude are  */
    specdat->latitude  = {lat:.4f};  /*   in DECIMAL DEGREES, not Deg/Min/Sec */
    specdat->timezone  =  -5.0;   /* Eastern time zone, even though longitude would
                                                suggest Central.  We use what they use.
                                                DO NOT ADJUST FOR DAYLIGHT SAVINGS TIME. */

    /*  */
    
    specdat->year      = {year:d};
    specdat->month     = {month:d};
    specdat->day       = {day:d};

    /*  */

    specdat->hour      = {hour:d};
    specdat->minute    = {minute:d};
    specdat->second    = {second:d};

    /* Let's assume that the temperature is 27 degrees C and that
       the pressure is 1006 millibars.  The temperature is used for the
       atmospheric refraction correction, and the pressure is used for the
       refraction correction and the pressure-corrected airmass. */

    specdat->temp      = {temp:.2f};
    specdat->press     = {press:.2f};

    /* We will use the first set of units (energy vs wavelength) */
    
    specdat->units     = 1;
    
    /* Tau at 500 nm will be assumed to be 0.2, and we will assume 1.36 cm of
       atmospheric water vapor in a vertical column. */
       
    specdat->tau500    = 0.2;
    specdat->watvap    = 1.36;
    specdat->ozone     = -1.0;  /* -1.0 input forces an internal calculation */

    /* Finally, we will assume that you have a flat surface facing southeast,
       tilted at latitude. */

    specdat->tilt      = specdat->latitude;  /* Tilted at latitude */
    specdat->aspect    = 180.0;       /* 0 deg. = facing north; 180 deg. = facing south */

    /* run the model */
    S_spectral2 ( specdat );

    /* after modifying the model code to save SZA, print it out here too */
    /*printf ( "%.4f", specdat->zenref )*/

    /*printf ( "micron  direct  diffuse\n" );*/
    int i;
    for (i = 0; i < 122; i++) {{

        printf ( "%.3e,%.6e,%.6e\n", 
            specdat->specx[i], specdat->specdir[i], specdat->specdif[i] );

    }}


        return ( 0 );
}}
"""


class spctral2():

    def __init__(self, 
        lat=45.0, lon=-80.0, 
        year=2014, month=7, day=1,
        hour=12, minute=0, second=0,
        temp=20, press=1000,
        casename='test', ID='001', 
        )

        self.lat = lat
        self.lon = lon

        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second
        
        self.temp = temp
        self.press = press

        # ozone, water, etc

        #> get ready for output
        self.casename = casename  # group of runs
        self.ID = ID  # identifier for one run of the group
        dirbases = ['img', 'raw', 'corrected']
        dirroot = './out/{:s}'.format(casename)
        self.outdir = dirroot
        for dirbase in dirbases:
            os.makedirs(dirroot + dirbase, exist_ok=True)

        self.raw_ofname = '{:s}/raw/{:s}.csv'.format(self.outdir, self.ID)
        self.corrected_ofname = '{:s}/corrected/{:s}.csv'.format(self.outdir, self.ID)
        


    def run(self)
        """ """
        
        #> create C code
        C_code = C_template.format(lon=lon, lat=lat,
                year=year, month=month, day=day,
                hour=hour, minute=minute, second=second,
                temp=Tair_C, press=p_mb)
        #print C_code

        cfname = 'tmp.c'
        with open(cfname, 'w') as f:
            f.writelines(C_code)

        #> compile C code

        xfname = 'tmp.out'
        cmd = 'icc {cfname:s} spectrl2.c solpos.c -o {xfname:s}'.format(cfname=cfname, xfname=xfname) 
        call(cmd.split())

        #> run the C code
        #  catching output in ofname

        #cmd2 = './{xfname:s} > {ofname:s}'.format(xfname=xfname, ofname=ofname)
        cmd2 = './{xfname:s}'.format(xfname=xfname)
        with open(self.raw_ofname, 'w') as f:
            call(cmd2, stdout=f)

        

    def correct(self, 
        measured_val=None, measured_bnds=()):
        """ 
        Check for negative values and report ToD and PPFD

        Change negative values (unphysical)
        and unrealistically high values to 0
        as well as nan's

        Correct total irradiance using a broadband measurement
          PAR, total, Vis, UV, etc
          Bounds of the measured irradiance value in microns!

        """

        #files = glob('{:s}/raw/*.csv'.format(self.outdir))


        wl, Idr, Idf = np.loadtxt(self.raw_ofname, delimiter=',', unpack=True)

        Idr[np.isnan(Idr)] = 0.0
        Idf[np.isnan(Idf)] = 0.0

        # correct for unrealistic and unphysical values
        Idr[(Idr < 0) & (Idr > 2000)] = 0
        Idf[(Idf < 0) & (Idf > 2000)] = 0  # < pretty sure only diffuse have negative values in the SPCTRAL2 output


        dwl = np.diff(wl)
        dwl = np.append(dwl, dwl[-1])  # or np.nan

        wl_meas = (wl >= measured_bnds[0]) & (wl <= measured_bnds[1])

        I_meas_sp2_all = (Idr+Idf)[wl_meas]  # all bands within the measured broad band
        dwl_meas = dwl[wl_meas]

        I_meas_sp2 = ( I_meas_sp2_all*dwl_meas / (6.626e-34*2.998e8/(wl[PAR_wl]*1e-6)) ).sum() / 6.022e23 * 1e6

        cf = measured_val / I_meas_sp2  # correction factor using measured broadband irradiance value

        Idr2 = Idr * cf
        Idf2 = Idf * cf 

        #> recheck the broadband calculation with sp2 data to make sure correction worked

        I_meas_sp2_all_2 = (Idr2+Idf2)[wl_meas]

        I_meas_sp2_2 = ( I_meas_sp2_all_2*dwl_meas / (6.626e-34*2.998e8/(wl[PAR_wl]*1e-6)) ).sum() / 6.022e23 * 1e6

        assert( np.all(np.isclose(PAR_spectral2_2, PAR_measured)) )


        #> save corrected spectra
        
        fnew = self.corrected_ofname

        a = np.vstack((wl, dwl, Idr2, Idf2)).T
        print a.shape
        np.savetxt(fnew, a, delimiter=',',
            fmt=['%.3e', '%.3e', '%.6e', '%.6e'])












