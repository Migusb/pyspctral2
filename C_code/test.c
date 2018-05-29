
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "spectrl2.h"   /* <-- include for spectrl2.  */


int main ( )
{
	struct specdattype spdat, *specdat;  /* spectral2 data structure */
	
	/* point to the specral2 structure */
	specdat = &spdat;
	
	/* Initialize the data structure -- you can omit if ALL structure vaiables 
	   		are assigned elswhere. Defaults in initialization can be overridden by 
	   		reassigning below. */
	   
	S_spec_init (specdat);
    
    /*  */
    
    specdat->longitude = -50.0000;  /* Note that latitude and longitude are  */
    specdat->latitude  = 40.0000;  /*   in DECIMAL DEGREES, not Deg/Min/Sec */
    specdat->timezone  =  -5.0;   /* Eastern time zone, even though longitude would
                                		suggest Central.  We use what they use.
                                		DO NOT ADJUST FOR DAYLIGHT SAVINGS TIME. */

    /*  */
    
    specdat->year      = 2014;
    specdat->month     = 7;
    specdat->day       = 14;

    /*  */

    specdat->hour      = 12;
    specdat->minute    = 0;
    specdat->second    = 0;

    /* Let's assume that the temperature is 27 degrees C and that
       the pressure is 1006 millibars.  The temperature is used for the
       atmospheric refraction correction, and the pressure is used for the
       refraction correction and the pressure-corrected airmass. */

    specdat->temp      = 27.00;
    specdat->press     = 1000.00;

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
    specdat->aspect    = 0.0;       /* 0 deg. = facing north */

    /* run the model */
    S_spectral2 ( specdat );

    /*printf ( "micron  direct  diffuse\n" );*/
    int i;
    for (i = 0; i < 122; i++) {

        printf ( "%.3e,%.6e,%.6e\n", 
            specdat->specx[i], specdat->specdir[i], specdat->specdif[i] );

    }


	return ( 0 );
}
