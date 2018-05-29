
"""

run the C version of SPCTRAL2 for each of the env data time steps 
for Chestnut Ridge 5-7 July sim


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

    /*printf ( "micron  direct  diffuse\n" );*/
    int i;
    for (i = 0; i < 122; i++) {{

        printf ( "%.3e,%.6e,%.6e\n", 
            specdat->specx[i], specdat->specdir[i], specdat->specdif[i] );

    }}


	return ( 0 );
}}
"""

#> location of Chesnut Ridge canopy site
lon = -84.2875
lat = 35.9583

#year = 2014
#month = 7
#day = 14
#
#day_first = 5
#num_days = 3
#
#dt0 = dt.datetime(year, month, day_first, 0, 0, 0) 
#
#dts = [dt0 + dt.timedelta(minutes=x) for x in range(30, (num_days*24*2+1)*30, 30)]
#
#for i, dti in enumerate(dts):
#    print i, dti


#> load time info, pressure, air T, PPFD from the env data file

env_data_f = '/storage/home/zlm1/work/access/zm/data/CH14_186-188.dat'

#sdate, stime, Tair_C, p_mb, PPFD = np.genfromtxt(env_data_f, 
env_data = np.genfromtxt(env_data_f, 
    skip_header=14, skip_footer=0, 
    usecols=(0, 1, 2, 4, 7),
    #names=['date', 'time', 'tac', 'pmb', 'ppfd'],
    dtype='S10,S8,f8,f8,f8',
    #unpack=True,
    )
#^ should just do 2 loads. this structure is not as easy to use as actual np arrays

#cmd0 = 'rm ./CR_201407_3day/raw/*.csv'

def run_spctral2_all():

    for i, line in enumerate(env_data):
        
        sdate, stime, Tair_C, p_mb, PPFD = line  # unpack

        sdt = sdate + ' ' + stime
        dtline = dt.datetime.strptime(sdt, '%Y-%m-%d %H:%M:%S')
        #print dtline

        year   = dtline.year
        month  = dtline.month
        day    = dtline.day
        hour   = dtline.hour
        minute = dtline.minute
        second = dtline.second

        #print year, month, day, hour, minute, second
     

        C_code = C_template.format(lon=lon, lat=lat,
                year=year, month=month, day=day,
                hour=hour, minute=minute, second=second,
                temp=Tair_C, press=p_mb)

        #print C_code

        cfname = 'tmp.c'
        with open(cfname, 'w') as f:
            f.writelines(C_code)

        xfname = 'tmp.out'
        cmd = 'icc {cfname:s} spectrl2.c solpos.c -o {xfname:s}'.format(cfname=cfname, xfname=xfname)  
        call(cmd.split())  

        ofname = './CR_201407_3day/raw/CRspectra{:03d}.csv'.format(i)  # catch output in this file
        #cmd2 = './{xfname:s} > {ofname:s}'.format(xfname=xfname, ofname=ofname)
        cmd2 = './{xfname:s}'.format(xfname=xfname)
        with open(ofname, 'w') as f:
            call(cmd2, stdout=f)


#> check for negative values and report tod and PPFD 
#  change negative values to 0 ? or just make them positive?

def correct_values_all():

    files = glob('./CR_201407_3day/raw/*.csv')

    for i, f in enumerate(sorted(files)):

        wl, Idr, Idf = np.loadtxt(f, delimiter=',', unpack=True)

        #drl0 = Idr < 0
        #dfl0 = Idf < 0
        #if np.any(drl0) or np.any(dfl0):
        #    print 'neg. values found. i = {:d}'.format(i)
        #    print '  stime = {:s}, PPFD = {:.1e}'.format(env_data[i][1], env_data[i][4])
        #    print '    direct: {:d} vals lt 0'.format(Idr[drl0].size)
        #    #print ','.join(np.char.mod('%.2e', direct[drl0]))
        #    print '    diffuse: {:d} vals lt 0'.format(Idf[dfl0].size)
        #    #print ','.join(np.char.mod('%.2e', diffuse[dfl0]))


        #> compute raw SPCTRAL2 PAR in umol photons / m^2 / s
        #  to compare with the measurement and create correction factor

        Idr[np.isnan(Idr)] = 0.0
        Idf[np.isnan(Idf)] = 0.0

        # correct for unrealistic and unphysical values
        Idr[(Idr < 0) & (Idr > 2000)] = 0
        Idf[(Idf < 0) & (Idf > 2000)] = 0  # < pretty sure only diffuse have negative values in the SPCTRAL2 output

        dwl = np.diff(wl)
        dwl = np.append(dwl, dwl[-1])  # or np.nan
        PAR_wl = (wl >= 0.4) & (wl <= 0.7)

        I_PAR = (Idr+Idf)[PAR_wl]  # irradiance dist (W/m^2/um)
        dwl_PAR = dwl[PAR_wl]  # um
        
        PAR_spectral2 = ( I_PAR*dwl_PAR / (6.626e-34*2.998e8/(wl[PAR_wl]*1e-6)) ).sum() / 6.022e23 * 1e6
        PAR_measured = env_data[i][4]

        cf = PAR_measured / PAR_spectral2  # correction factor, to be applied across all bands

        #print 'i = {:d}'.format(i)
        #print '  PAR {:.2e} vs {:.2e}'.format(PAR_spectral2, PAR_measured)
        

        Idr2 = Idr * cf  # applying correction factor
        Idf2 = Idf * cf

        #> check PAR calculation again to make sure correction worked

        I_PAR_2 = (Idr2+Idf2)[PAR_wl]

        PAR_spectral2_2 = ( I_PAR_2*dwl_PAR / (6.626e-34*2.998e8/(wl[PAR_wl]*1e-6)) ).sum() / 6.022e23 * 1e6

        assert( np.all(np.isclose(PAR_spectral2_2, PAR_measured)) )


        #> save corrected spectra
        fbn = os.path.basename(f)
        fnew = './CR_201407_3day/corrected/{:s}'.format(fbn)
        
        a = np.vstack((wl, dwl, Idr2, Idf2)).T
        print a.shape
        np.savetxt(fnew, a, delimiter=',', 
            fmt=['%.3e', '%.3e', '%.6e', '%.6e'])    



if __name__ == '__main__':

    #run_spctral2_all()

    correct_values_all()


