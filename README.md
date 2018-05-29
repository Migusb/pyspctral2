# pyspctral2
Basic Python wrapper for functionality of the SPCTRAL2 atmospheric radiative transfer C code

## SPCTRAL2

The Bird Simple Spectral Model v2 ([SPCTRAL2](http://rredc.nrel.gov/solar/models/spectral/)) is a simple atmospheric radiative transfer model that predicts spectral irradiance at the surface.
* 122 bands, in region 0.3&ndash;4.0 &mu;m, with higher spectral resolution in UV and Vis
* Predicts diffuse light and direct solar beam
* Processes modeled using semi-empirical exponential extinction:
  * Ozone
  * Water vapor
  * Refraction

### Inputs

Required
* Location (lat, lon)
* Date, time, UTC offset
* Angle (tilt and aspect) of flat collector


