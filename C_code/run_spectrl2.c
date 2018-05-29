/*============================================================================
*
*    NAME:  spectest.c
*
*    PURPOSE:  Exercises NREL's Simple Spectral model, 'spectrl2.c'.
*
*        Based on the SERI (now NREL) technical report SERI/TR-215-2436,
*        "Simple Solar Spectral Model for Direct and Diffuse Irradiance
*        on Horizontal and Tilted Planes at the Earth's Surface for 
*        Cloudless Atmospheres", by
*            R. Bird & C. Riordan
*
*    'spectrl2.c' contains:
*
*        int S_spectral (void)
*
*        INPUTS:
*            alpha      DEFAULT 1.14, power on Angstrom turbidity
*            aspect     DEFAULT 180.0 (South; N=0, E=90, W=270),
*                            aspect (azimuth angle) of collector panel.
*            assym      DEFAULT 0.65, aerosol assymetry factor
*			 day		Day of month (May 27 = 27, etc.) 
*			 hour		Hour of day, 0 - 23
*			 minute		Minute of hour, 0 - 59
*			 month		Month number (Jan = 1, Feb = 2, etc.)
*            latitude	Latitude, degrees north (south negative)
*            longitude	Longitude, degrees east (west negative)
*            ozone      DEFAULT -1.0 will force an internal calculation;
*                            ozone amount (atmospheric cm)
*            press		Surface pressure, millibars
*			 second		Second of minute, 0 - 59
*            spcrfl[6]  DEFAULT 0.2 (all), ground reflectivities
*            spcwvr[6]  DEFAULT { 0.3, 0.7, 0.8, 1.3, 2.5, 4.0 }
*                            wavelength regions (microns) for reflectivities
*            tau500     DEFAULT -1.0 (missing);
*                            aerosol optical depth at 0.5 microns, base e
*            tilt       DEFAULT 0.0, tilt angle of collector panel;
*                            if tilt > 180.0, sun-tracking is assumed.
*            timezone	Time zone, east (west negative)
*            temp		Ambient dry-bulb temperature, degrees C
*            units      Output units:
*                            1 = irradiance (W/sq m/micron) 
*                                  per wavelength (microns)
*                            2 = photon flux (10.0E+16 /sq cm/s/micron) 
*                                  per wavelength (microns)
*                            3 = photon flux density (10.0E+16 /sq cm/s/eV)
*                                  per energy (eV)
*            watvap     DEFAULT -1.0 (missing);
*                            precipitable water vapor (cm)
* 			 year		4-digit year
* 
*            specdif[122]    Diffuse spectrum on panel (see units)
*            specdir[122]    Direct normal spectrum (see units)
*            specetr[122]    Extraterrestrial spectrum (W/sq m/micron)
*            specglo[122]    Global spectrum on panel (see units)
*            specx[122]      X-coordinate (wavelength or energy; see units)
*
*    Usage:
*         In calling program, just after other 'includes', insert:
*
*              #include "spectrl2.h"
*
*         "solpos.c" must be linked as well as "spectrl2.c", and its inputs
*         must be satisfied as well (see its documentation).
*
*
*    Martin Rymes
*    National Renewable Energy Laboratory
*    25 March 1998
*----------------------------------------------------------------------------
*
* Modifications:
* -------------
*
*	18 March 2004: Variables now passed to functions through structure.
*
*		Steve Wilcox
*		Mary Anderberg
*		National Renewable Energy Laboratory
*
*----------------------------------------------------------------------------*/
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
    
    /* I use Atlanta, GA for this example */
    
    specdat->longitude = -84.43;  /* Note that latitude and longitude are  */
    specdat->latitude  =  33.65;  /*   in DECIMAL DEGREES, not Deg/Min/Sec */
    specdat->timezone  =  -5.0;   /* Eastern time zone, even though longitude would
                                		suggest Central.  We use what they use.
                                		DO NOT ADJUST FOR DAYLIGHT SAVINGS TIME. */

    /* The date is 22-July-2003 */
    
	specdat->year      = 1999;
    specdat->month     =    7;
    specdat->day       =   22;

    /* The time of day (STANDARD time) is 9:45:37 */

    specdat->hour      = 9;
    specdat->minute    = 45;
    specdat->second    = 37;

    /* Let's assume that the temperature is 27 degrees C and that
       the pressure is 1006 millibars.  The temperature is used for the
       atmospheric refraction correction, and the pressure is used for the
       refraction correction and the pressure-corrected airmass. */

    specdat->temp      =   27.0;
    specdat->press     = 1006.0;

    /* We will use the first set of units (energy vs wavelength) */
    
    specdat->units     = 1;
    
    /* Tau at 500 nm will be assumed to be 0.2, and we will assume 1.36 cm of
       atmospheric water vapor in a vertical column. */
       
    specdat->tau500    = 0.2;
    specdat->watvap    = 1.36;

    /* Finally, we will assume that you have a flat surface facing southeast,
       tilted at latitude. */

    specdat->tilt      = specdat->latitude;  /* Tilted at latitude */
    specdat->aspect    = 0.0;       /* 0 deg. = facing north */

    printf ( "\n" );
    printf ( "***** TEST S_spectral2: *****\n" );
    printf ( "\n" );
    
	S_spectral2 ( specdat );

    /*printf ("\n" );
    printf ( "micron  (NREL) ETRglo (NREL) global (NREL) direct (NREL) diffuse (NREL)\n" );
    printf ( "micron  (NREL) ETRglo (NREL) global (NREL) direct (NREL) diffuse (NREL)\n" );
    printf ( "%6.3f  0.300 %6.1f  518.7 %6.1f    2.7 %6.1f    0.9 %6.1f    1.8\n",
        specdat->specx[  0], specdat->specetr[  0], specdat->specglo[  0], specdat->specdir[  0], specdat->specdif[  0] );
    */

    printf ( "micron  direct  diffuse\n" );
    int i;
    for (i = 0; i < 122; i++) {

        printf ( "%.3e,%.6e,%.6e\n", 
            specdat->specx[i], specdat->specdir[i], specdat->specdif[i] );

    }


	return ( 0 );
}
