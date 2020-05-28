
"""
Python module to run the SPCTRAL2 C code and catch output


typical call sequence:
.. code-block::
   import pyspctral2 as psp2

"""

import datetime as dt
from glob import glob
import os
#from pathlib2 import Path  # for py2 compatability
from subprocess import call

import numpy as np
import xarray as xr


_this_dir = os.path.dirname(os.path.realpath(__file__))
C_code_path = '{:s}/C_code'.format(_this_dir)  # loc of SPCTRAL2 and solpos C code
C_template_path = '{:s}/run_spctral2_template.c'.format(C_code_path)
cwd = os.getcwd()

micron = u'\u03BCm'

with open(C_template_path, 'r') as f:
    C_template = f.read()



SPCTRAL2_MODEL_REQUIRED_INPUT_PARAMETER_NAMES = [
    "lat", "lon", 
    "year", "month", "day", "hour", "minute", "second",
    "utcoffset",
    "temp", "press",
    "tau500", "watvap", "ozone",
    "tilt", "aspect",
]


# just to see if it will run or not... (not for correctness)
def test_run_csp2_py():
    from C_code_cython.csp2_py import run_spctral2

    return run_spctral2(lat=40.0, lon=-105.0, 
        year=2014, month=7, day=1,
        hour=12, minute=0, second=0,
        utcoffset=-5.,
        temp=27, press=1010,
        tau500=0.27, watvap=1.42, ozone=-1.0)


def _is_tool(name):
    """Check whether `name` is on PATH and marked as executable.
    
    ref: https://stackoverflow.com/a/34177358
    """

    from shutil import which

    return which(name) is not None


def _recompile_run(m):
    """
    Use a modified version of the runner program 
    provided with the SPCTRAL2 C distribution (`spectest.c`).

    This requires modifying the C code of the runner program and recompiling it 
    every time. The runner program writes the data to a text file.

    PARAMETERS
    ----------
    m : model object

    RETURNS
    -------
    None
        since for now we pass and modify the model object

    """

    # hack for now
    self = m

    #> create C code
    C_code = C_template.format(lon=self.lon, lat=self.lat,
        year=self.year, month=self.month, day=self.day,
        hour=self.hour, minute=self.minute, second=self.second,
        utcoffset=self.utcoffset,
        temp=self.temp, press=self.press, 
        tau500=self.tau500, watvap=self.watvap, ozone=self.ozone,
        )
    #print C_code

    os.chdir(C_code_path)

    cfname = 'tmp.c'
    with open(cfname, 'w') as f:
        f.writelines(C_code)

    #> compile C code
    # TODO: allow compiler choice as fn param
    
    # for now use gcc but check
    if not _is_tool("gcc"):
        raise Exception(
            "GNU C compiler `gcc` not found. "
            "Must install to use this SPCTRAL2 run method."
        )
    # for more general case should use `try: ...` and handle `OSError` if it arises
    # ref: https://stackoverflow.com/a/11210185
    
    xfname = 'tmp.out'
    #cmd = 'icc {cfname:s} spectrl2.c spectrl2.h solpos.c solpos00.h -o {xfname:s}'.format(cfname=cfname, xfname=xfname) 
    #cmd = 'gcc {cfname:s} spectrl2.c solpos.c -o {xfname:s}'.format(cfname=cfname, xfname=xfname) 
    cmd = 'gcc {cfname:s} spectrl2.c solpos.c -o {xfname:s} -lm'.format(cfname=cfname, xfname=xfname) 
    # ^ I guess I was using a quite old gcc before (linking math was not required)
    # https://stackoverflow.com/a/5005419

    call(cmd.split())
    os.chdir(cwd)

    #> run the C code
    #  catching output in ofname

    #cmd2 = './{xfname:s} > {ofname:s}'.format(xfname=xfname, ofname=ofname)
    cmd2 = '{c:s}/{xfname:s}'.format(xfname=xfname, c=C_code_path)
    with open(self.raw_ofname, 'w') as f:
        call(cmd2, stdout=f)

    #> load raw data from the SPCTRAL2 run 
    #  first row in raw file is NREL SPA computed solar zenith angle
    #  the rest of the rows are three cols: wavelength, direct, diffuse
    with open(self.raw_ofname, 'r') as f:
        sza = float(f.readline().rstrip())
    wl, Idr, Idf = np.loadtxt(self.raw_ofname, delimiter=',', skiprows=1, unpack=True)

    return sza, wl, Idr, Idf



class model():
    """Wrapper for the SPCTRAL2 model (C version).
    
    Note that there are no runtime options in the original model,
    so we compile every time (not very efficient...)
    but luckily it takes very little time to run. 

    Initialization parameters
    -----------------------------
    lat, lon : float
    year, month, day, hour, minute, second : int
        local time!
    utcoffset : int
        UTC offset (negative for west)
    tau500 : float
        aerosol optical depth at 0.5 um (500 nm), base e
        set to 0.0 to ignore (-1.0 to silently fail with warning in output file...)
        default value (from original C example run script) = 0.2, from Excel version = 0.27 (midrange)
        typical range [0.05, 0.55] for clear sky
    watvap : float
        column precipitable water vapor (cm)
        set to 0.0 to ignore (-1.0 to silently fail with warning in output file...)
        default value (from original C example run script) = 1.36, from Excel version = 1.42
        typical range [0.3, 6]
    ozone : float
        total column ozone (cm)
        set to -1.0 for an internal (SPCTRAL2) estimation based on lat/lon and time of day/year (Heuklon 1978)
        default value (from original C example run script) = -1.0
    casename : str
        name of the case, used to create a directory for saving the outputs
        default = 'test'
    ID : str
        identifier used for the saved output filename
        default = '001'

    """

    def __init__(self, 
        lat=40.0, lon=-105.0, 
        year=2014, month=7, day=1,
        hour=12, minute=0, second=0,
        utcoffset=-5.,
        temp=27, press=1010,
        tau500=0.27, watvap=1.42, ozone=-1.0, 
        casename='test', ID='001',
        sza_check=None,
        cc_correct=False, cc_p1=1.0, cc_p2=0.3,
        ):

        self.lat = lat
        self.lon = lon

        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second
        self.utcoffset = utcoffset
        self.dt = dt.datetime(year, month, day, hour, minute, second)  # could add timezone info?
 
        self.temp = temp    # near-sfc air temp (deg. C)
        self.press = press  # near-sfc air pressure (mb)

        self.cc_correct = cc_correct
        self.cc_p1 = cc_p1
        self.cc_p2 = cc_p2

        # ozone, water, etc
        self.tau500 = tau500
        if self.cc_correct:
            self.watvap = watvap * 1.2 * 1.2 * self.cc_p1
            #^ climatologically higher watvap for cloudy conditions + 
            #  increase in water vapor path due to multiple scattering within cloud
        else:
            self.watvap = watvap
        self.ozone = ozone

        self.sza_check = sza_check

        #> get ready for output

        self.casename = casename  # group of runs
        self.ID = ID  # identifier for one run of the group
        dirbases = ['img', 'raw', 'corrected']
        dirroot = './out/{:s}'.format(casename)
        
        # TODO: check if desired ID already exists in the out dir and give warning or error
        # TODO: instead of 'out', more descriptive name (like pyspctral2_out) 
        #       for output directory (needs to be changed multiple places)

        self.outdir = dirroot

        for dirbase in dirbases:
            os.makedirs('{:s}/{:s}'.format(dirroot, dirbase), exist_ok=True)  # exist_ok only in py3 version
            #ps = '{:s}/{:s}'.format(dirroot, dirbase)  # py2 exist_ok workaround using pathlib2
            #Path(ps).mkdir(parents=True, exist_ok=True)

        self.raw_ofname = '{:s}/raw/{:s}.csv'.format(self.outdir, self.ID)
        self.corrected_ofname = '{:s}/corrected/{:s}.csv'.format(self.outdir, self.ID)
        


    def run(self):
        """ 
        Run the SPCTRAL2 C program using the selected method (C code or Cython). 
        """
        
        # run using selected runner
        sza, wl, Idr, Idf = _recompile_run(self)  # only this one is set up


        #> check sza (if check value provided)
        #  first row in raw file is solar zenith angle
        if self.sza_check:
            # with open(self.raw_ofname, 'r') as f:
            #     sza_sp2 = float(f.readline().rstrip())
            sza_sp2 = sza
            if abs(sza_sp2 - self.sza_check) > 0.1:  # should be very close
                print('SZA issue:')
                print('  sp2  : {:.3f}'.format(sza_sp2))
                print('  check: {:.3f}'.format(self.sza_check))
            else:
                print('SZA matches check.')

        

    def correct(self, 
        measured_units='E', measured_val=None, measured_bnds=(0.4, 0.7),
        total_solar=None,
        plot=False, save_plot=False,
        ):
        """Apply corrections to raw SPCTRAL2 outputs

        1. Check for negative values and report ToD and PPFD

        2. Change negative values (unphysical) and unrealistically high values to 0
           as well as nan's

        3. Correct total irradiance using a total (direct+diffuse) broadband measurement
           e.g., PAR, total, Vis, UV, etc
           Could do better job if we had separate direct and diffuse measurements...

        Parameters
        ----------
        measured_units : str
            'E' (W/m^2) or 'photons' (photon flux density)
        measured_val : float
            A measured broadband irradiance value (umol photons / m^2 / s)
        measured_bnds : tuple
            Wavelength band bounds for the measured irradiance value (um) 
        total_solar : float
            A measurement of total solar irradiance (direct+diffuse; W/m^2)
            We can use this to correct the regions outside the measured band (e.g., PAR)
            Only use this in addition to a measured band value!


        """

        #files = glob('{:s}/raw/*.csv'.format(self.outdir))

        #> load raw data from the SPCTRAL2 run 
        #  first row in raw file is NREL SPA computed solar zenith angle
#        with open(self.raw_ofname, 'r') as f:
#            sza = f.readline().rstrip()
        wl, Idr, Idf = np.loadtxt(self.raw_ofname, delimiter=',', skiprows=1, unpack=True)

        #> NaN -> 0
        Idr[np.isnan(Idr)] = 0.0
        Idf[np.isnan(Idf)] = 0.0

        #> correct for unrealistic and unphysical values
        Idr[(Idr < 0) | (Idr > 2000)] = 0
        Idf[(Idf < 0) | (Idf > 2000)] = 0  # < pretty sure only diffuse have negative values in the SPCTRAL2 output

        #drl0 = Idr < 0
        #dfl0 = Idf < 0
        #if np.any(drl0) or np.any(dfl0):
        #    print 'neg. values found. i = {:d}'.format(i)
        #    print '  stime = {:s}, PPFD = {:.1e}'.format(env_data[i][1], env_data[i][4])
        #    print '    direct: {:d} vals lt 0'.format(Idr[drl0].size)
        #    #print ','.join(np.char.mod('%.2e', direct[drl0]))
        #    print '    diffuse: {:d} vals lt 0'.format(Idf[dfl0].size)
        #    #print ','.join(np.char.mod('%.2e', diffuse[dfl0]))

        #> calculate band widths
        dwl = np.diff(wl)
        dwl = np.append(dwl, dwl[-1])  # or np.nan

        # TODO: compare diffuse/direct ratio with Spitters method for PAR (at least)

        # TODO: cloudiness correction (for low-level clouds, non-precip):
        # - 2 spectral modifications from Bird et al. (1997)
        # - also adjust watvap, 
        #   to include effect of increased path length due to multiple scattering in cloud
        # maybe 2 parameters, 
        #   for cloud thickness and cloud fraction in the area
        #   don't worry about cloud/hydrometeor types

        #> Bird et al. (1987) suggested cloud correction to diffuse spectrum
        #  but adding a parameter to increase/decrease the magnitude of the changes
        if self.cc_correct:
            wl_cc_1 = wl <= 0.55
            wl_cc_2 = (wl >= 0.50) & (wl <= 0.926)
            cc_1 = (wl[wl_cc_1] + 0.45)**(-1.0)  # first correction to Idf
            cc_2 = (1 + 0.07*3)                  # second correction
            
            # Bird et al. (1987) correction as is (only diffuse increased)
            Idf[wl_cc_1] *= cc_1 
            Idf[wl_cc_2] *= cc_2

            # attempting an E-conserving version
            # Idf[wl_cc_1] += Idr[wl_cc_1] * (1-1/cc_1) * self.cc_p1 #* 0.07/0.333
            # Idf[wl_cc_2] += Idr[wl_cc_2] * (1-1/cc_2) * self.cc_p1
            # #Idr[wl_cc_1] *= 1/cc_1 
            # #Idr[wl_cc_2] *= 1/cc_2
            # Idr[wl_cc_1] -= Idr[wl_cc_1] * (1-1/cc_1) * self.cc_p1
            # Idr[wl_cc_2] -= Idr[wl_cc_2] * (1-1/cc_2) * self.cc_p1

            # convert a fraction of direct to diffuse
            # based on input fraction related to how much the Sun is obscured 
            # i.e., transmissivity of the direct beam through cloud layer
            Idf += Idr * (1-self.cc_p2)
            Idr -= Idr * (1-self.cc_p2)
            #^ for now at all wavelengths equally


        #--------------------------------------------------------------------------------------------
        #> normalize region(s) to match measured values 

        assert( len(measured_bnds) == 2 and measured_bnds[0] < measured_bnds[1] )
        wl_meas = (wl >= measured_bnds[0]) & (wl <= measured_bnds[1])
        wl_meas_not = ~wl_meas

        I_meas_sp2_all = (Idr+Idf)[wl_meas]  # spectral irradiances for all sub-bands within the measured broadband
        dwl_meas = dwl[wl_meas]              # corresponding bandwidths

        I_meas_not_sp2_all = (Idr+Idf)[wl_meas_not]
        dwl_meas_not = dwl[wl_meas_not]

        if measured_units == 'photons':
            # convert SPCTRAL2 sub-band irradiances in W/m^2/um to broadband in photon flux units (photons/m^2/s)
            I_meas_sp2 = ( I_meas_sp2_all*dwl_meas / (6.626e-34*2.998e8/(wl[wl_meas]*1e-6)) ).sum() / 6.022e23 * 1e6
        elif measured_units == 'E':  # E units (W/m^2)
            I_meas_sp2 = ( I_meas_sp2_all*dwl_meas ).sum()
        else:
            pass  # should raise error
        I_meas_not_sp2 = ( I_meas_not_sp2_all*dwl_meas_not ).sum()

        if measured_val:
            cf_band = measured_val / I_meas_sp2  # correction factor using measured broadband irradiance value
            if total_solar:
                #cf = np.ones_like(wl)
                cf = np.ones(wl.shape)  # to please pylint
                cf[wl_meas] = cf_band
                cf[wl_meas_not] = (total_solar-measured_val) / I_meas_not_sp2
            else:
                cf = cf_band
        else:
            cf = 1.0

        Idr2 = Idr * cf
        Idf2 = Idf * cf 

        #> recheck the broadband calculation with sp2 data to make sure correction worked
        I_meas_sp2_all_2 = (Idr2+Idf2)[wl_meas]
        if measured_units == 'photons':
            I_meas_sp2_2 = ( I_meas_sp2_all_2*dwl_meas / (6.626e-34*2.998e8/(wl[wl_meas]*1e-6)) ).sum() / 6.022e23 * 1e6
        elif measured_units == 'E':
            I_meas_sp2_2 = ( I_meas_sp2_all_2*dwl_meas ).sum()
        else:
            pass  # should raise error

        if measured_val:
            assert( np.all(np.isclose(I_meas_sp2_2, measured_val)) )

        #> check total solar?
        if total_solar:
            total_solar_check = ( (Idr2+Idf2)*dwl ).sum()
            assert( np.all(np.isclose(total_solar, total_solar_check)) )
            print('total solar checks out.')
        
        #--------------------------------------------------------------------------------------------

        #> save corrected spectra
        fnew = self.corrected_ofname
        a = np.vstack((wl, dwl, Idr2, Idf2)).T
        #print a.shape
        np.savetxt(fnew, a, delimiter=',',
            fmt=['%.3e', '%.3e', '%.6e', '%.6e'])

        #> plot if desired
        if plot == True:
            plot_spectrum(casename=self.casename, ID=self.ID, output_type='corrected',
                title='', title_left='', title_right='',
                save=save_plot)


            
def combine(casename='blah', time=[], info=''):
    """
    For given casename, combine all of the individual output files and data
    into one netCDF file

    time: assumed list of datetime objs
        used to timestamp the spectra

    Ultimately would like to have 
      coords: wl, time
      vars: dwl, Idr, Idf, sdatetime, T, p, ozone, watvap, tau500
    """

    #> get lists of files and sort
    files_raw = glob('./out/{:s}/raw/*.csv'.format(casename))
    files_corrected = glob('./out/{:s}/corrected/*.csv'.format(casename))
    files_raw.sort()
    files_corrected.sort()

    assert( len(files_raw) == len(time) and len(files_raw) == len(files_corrected) )
    nt = len(time)
    nwl = 122

    #> time variables
    sdt = [dtx.strftime('%Y-%m-%d %H:%M:%S') for dtx in time]

    #> get sza from 1st line of the raw files
    #  SPCTRAL2 calculates it using the NREL SPA
    sza = []
    for i in range(len(files_raw)):
        with open(files_raw[i], 'r') as f:
            szai = float(f.readline().rstrip())
        sza.append(szai)

    #> build Idr and Idf arrays
    Idr = np.zeros((nt, nwl))
    #Idf = np.zeros_like(Idr)
    Idf = np.zeros(Idr.shape)  # to please pylint
    for i, f in enumerate(files_corrected):
        #print i, f
        wl, dwl, Idri, Idfi = np.loadtxt(f, delimiter=',', unpack=True)
        Idr[i,:] = Idri
        Idf[i,:] = Idfi    

    #> create dataset
    dset = xr.Dataset(\
        coords={'wl': ('wl', wl, {'units': micron, 'longname': 'wavelength'}), 
                'time': time}, 
        data_vars={\
            'Idr': (['time', 'wl'], Idr, {'units': 'W m^-2 um^-1', 'longname': 'direct beam solar irradiance'}),
            'Idf': (['time', 'wl'], Idf, {'units': 'W m^-2 um^-1', 'longname': 'diffuse solar irradiance'}),
            'dwl': ('wl', dwl, {'units': 'um', 'longname': 'wavelength band width'}),
            'sdt': ('time', sdt, {'longname': 'datetime string'}),
            'sza': ('time', sza, {'units': 'deg.', 'longname': 'solar zenith angle'}),
            },
        attrs={'casename': casename, 'info': info},
        )

    #> save dataset
    ofname = './out/{:s}.nc'.format(casename)
    dset.to_netcdf(ofname)




def plot_spectrum(casename='test', ID='001', output_type='raw',
    title='', title_left='', title_right='',
    save=False):
    """Plot spectrum from an output CSV file, given the filename.


    Returns
    -------
    pyplot figure handle
    """
    import matplotlib.pyplot as plt

    ofpath = './out/{casename:s}/{output_type:s}/{ID:s}.csv'.format(\
        casename=casename, output_type=output_type, ID=ID)
    
    if output_type == 'raw':
        wl, Idr, Idf = np.loadtxt(ofpath, delimiter=',', skiprows=1, unpack=True)
    elif output_type == 'corrected':
        wl, _, Idr, Idf = np.loadtxt(ofpath, delimiter=',', skiprows=1, unpack=True)
    else:
        print('invalid output_type selection')  # should raise error instead...
        return False

    # dwl = np.diff(wl)
    # dwl = np.append(dwl, dwl[-1])  # or np.nan

    f1, a = plt.subplots(figsize=(8, 4))

    a.plot(wl, Idr, label='direct')
    a.plot(wl, Idf, label='diffuse')
    a.plot(wl, Idr+Idf, label='total')
    
    a.set_xlim((0.3, 2.6))  # SPCTRAL2 limits

    a.set_title(title)
    a.set_title(title_left, loc='left')
    a.set_title(title_right, loc='right')
    a.set_xlabel('wavelength (μm)')
    a.set_ylabel('spectral irradiance (W m$^{-2}$ μm$^{-1}$)')
    
    a.legend()
    a.grid(True)
    f1.tight_layout()

    if save:
        outdir = './out/{casename:s}/img'.format(casename=casename)
        f1.savefig('{:s}/{:s}.png'.format(outdir, ID), dpi=150, 
            transparent=False, bbox_inches='tight', pad_inches=0.05)

    return f1
