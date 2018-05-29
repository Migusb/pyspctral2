
/* 
This is a C code template for use with Python string format

As such, curly braces are doubled,
in places where curly brace is intended within the C code


*/

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
 *                         are assigned elswhere. Defaults in initialization can be overridden by 
 *                                                 reassigning below. */
           
        S_spec_init (specdat);
    
    /*  */
    
    specdat->longitude = {lon:.4f};  /* Note that latitude and longitude are  */
    specdat->latitude  = {lat:.4f};  /*   in DECIMAL DEGREES, not Deg/Min/Sec */
    specdat->timezone  = {utcoffset: .1f};   /* -5: Eastern time zone, even though longitude would
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
 *        the pressure is 1006 millibars.  The temperature is used for the
 *               atmospheric refraction correction, and the pressure is used for the
 *                      refraction correction and the pressure-corrected airmass. */

    specdat->temp      = {temp:.2f};
    specdat->press     = {press:.2f};

    /* We will use the first set of units (energy vs wavelength) */
    
    specdat->units     = 1;
    
    /* Tau at 500 nm will be assumed to be 0.2, and we will assume 1.36 cm of
 *        atmospheric water vapor in a vertical column. */

    specdat->tau500    = 0.2;
    specdat->watvap    = 1.36;
    specdat->ozone     = -1.0;  /* -1.0 input forces an internal calculation */

    /* Finally, we will assume that you have a flat surface facing southeast,
 *        tilted at latitude. */

    specdat->tilt      = specdat->latitude;  /* Tilted at latitude */
    specdat->aspect    = 180.0;       /* 0 deg. = facing north; 180 deg. = facing south */

    /* run the model */
    S_spectral2 ( specdat );

    /* print solar zenith angle (calculated by solpos)
       units: not sure */
    printf ( "%.5f\n", specdat->zenref );

    /*printf ( "micron  direct  diffuse\n" );*/
    int i;
    for (i = 0; i < 122; i++) {{

        printf ( "%.3e,%.6e,%.6e\n", 
            specdat->specx[i], specdat->specdir[i], specdat->specdif[i] );

    }}


        return ( 0 );
}}


