/* 
* This is a C code template for use with Python string format method.
* 
* As such, curly braces are doubled,
* in places where curly brace is intended within the C code.
* 
* 
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
*         are assigned elswhere. Defaults in initialization can be overridden by 
*         reassigning below. */
   S_spec_init (specdat);
    
    /*  */
    specdat->longitude = {lon:.4f};  /* latitude and longitude in DECIMAL DEGREES, not Deg/Min/Sec (duh) */
    specdat->latitude  = {lat:.4f};
    specdat->timezone  = {utcoffset:.1f};   /* local time zone. DO NOT ADJUST FOR DAYLIGHT SAVINGS TIME.
*                                               e.g., -5: Eastern time zone */

    /*  */
    specdat->year      = {year:d};
    specdat->month     = {month:d};
    specdat->day       = {day:d};

    /*  */
    specdat->hour      = {hour:d};
    specdat->minute    = {minute:d};
    specdat->second    = {second:d};

    /* The temperature is used for the
 *     atmospheric refraction correction, and the pressure is used for the
 *     refraction correction and the pressure-corrected airmass. */
    specdat->temp      = {temp:.2f};
    specdat->press     = {press:.2f};

    /* We will use the first set of units (energy vs wavelength) 
*      this gives us spectral irradiances (W/m^2/um) and wavelengths (um) */
    specdat->units     = 1;
    
    /* In original example C run script: 
*        Tau at 500 nm assumed to be 0.2, and we will assume 1.36 cm of
*        atmospheric water vapor in a vertical column. */
    specdat->tau500    = {tau500:.3f};
    specdat->watvap    = {watvap:.3f};
    specdat->ozone     = {ozone:.3f};  /* -1.0 input forces an internal ozone estimate */

    /* Dectector plane specifications */
    specdat->tilt      = specdat->latitude;  /* Tilted at latitude */
    specdat->aspect    = 0.0;       /* 0 deg. = facing north; 180 deg. = facing south */
   /* For such a plane in southern hemisphere, tilt at latitude and facing north
*     should produce a plane parallel to the spherical Earth ground sfc at that latitude 
*     In the original provided run script, facing southeast (270 deg.) was used */ 

    /* Run the model */
    S_spectral2 ( specdat );

    /* Print solar zenith angle (calculated by solpos)
       units: not sure (appears to be deg.!) */
    printf ( "%.5f\n", specdat->zenref );

    /*printf ( "wavlength (micron)   direct (W/m^2/micron)  diffuse (W/m^2/micron)\n" );*/
    int i;
    for (i = 0; i < 122; i++) {{

        printf ( "%.3e,%.6e,%.6e\n", 
            specdat->specx[i], specdat->specdir[i], specdat->specdif[i] );

    }}


   return ( 0 );
}}


