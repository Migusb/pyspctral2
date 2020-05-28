
"""
Run the SPCTRAL2 C code.


typical call sequence:
.. code-block::
   import pyspctral2 as psp2

"""

import datetime
# from glob import glob
import os
import subprocess

import numpy as np
import xarray as xr


_THIS_DIR = os.path.dirname(os.path.realpath(__file__))
C_CODE_PATH = f'{_THIS_DIR:s}/C_code'  # loc of SPCTRAL2 and solpos C code
C_TEMPLATE_PATH = f'{C_CODE_PATH:s}/run_spctral2.c.tpl'
CWD = os.getcwd()
# TODO: upgrade to use pathlib (and check paths exist?)


SPCTRAL2_MODEL_REQUIRED_INPUT_PARAMETER_NAMES = [
    "lat", "lon", 
    "year", "month", "day", "hour", "minute", "second",
    "utcoffset",
    "temp", "press",
    "tau500", "watvap", "ozone",
    # "tilt", "aspect",  # < currently these values are hardcoded in the runners
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


def _within_range_check(value, bounds, *, mode="hard", name=""):
    # based on xarray-simlab's xsimlab.validators.in_bounds
    # see: https://github.com/benbovy/xarray-simlab/blob/master/xsimlab/validators.py
    # ^ could expand this one to work for np arrays too

    delim_left  = "["
    delim_right = "]"
    bounds_str = f"{delim_left}{bounds[0]}, {bounds[1]}{delim_right}"

    if (
        bounds[0] is not None
        and bounds[1] is not None
        and bounds[1] < bounds[0]
    ):
        raise ValueError(
            f"Invalid bounds {bounds_str}: "
            "Upper limit should be higher than lower limit."
        )

    out_lower = bounds[0] is not None and value <= bounds[0]

    out_upper = bounds[1] is not None and value >= bounds[1]

    if out_lower or out_upper:
        if name:
            msg_name = f" for variable `{name}`"
        else:
            msg_name = ""

        msg = f"Value {value} is out of bounds {bounds_str}{msg_name}"

        if mode == "hard":
            raise ValueError(msg)
        elif mode == "soft":
            print(f"Warning: {msg}")
        else:
            raise ValueError(f"`mode` {mode!r} is invalid.")



def validate_within_range_inclusive(x, hard_bounds=(None, None), soft_bounds=(), *, name=""):
    """Check if value `x` is within the bounds. 

    Bounds passed to `hard_bounds` raise an error; 
    those passed to `soft_bounds` a warning. 

    x : float
        value to check is within the bounds
    hard : 2-tuple, items can be float, int, or None
        raise ValueError if this one is not met
        default = unbounded (None, None)
    soft : 2-tuple, items can be float, int, or None; optional
        warning message instead of raising error
        should be a subset of the hard bounds (not extend outside)

    name : str, optional
        variable name to use in the outside-bounds message
    """
        
    # validate inputs
    def validate_bounds(bounds):
        NoneType = type(None)

        ok = (
            len(bounds) == 2 
            and all(isinstance(x, (float, int, NoneType)) for x in bounds)
        )
        if not ok:
            raise ValueError(
                f"{bounds!r} is an invalid bounds specification."
            )

    # check hard bounds
    validate_bounds(hard_bounds)
    _within_range_check(x, hard_bounds, mode="hard", name=name)

    # check soft bounds (optional)
    if soft_bounds:
        validate_bounds(soft_bounds)
        _within_range_check(x, soft_bounds, mode="soft", name=name)



def _recompile_run(sp2_inputs):
    """
    Use a modified version of the runner program 
    provided with the SPCTRAL2 C distribution (`spectest.c`).

    This requires modifying the C code of the runner program and recompiling it 
    every time. The runner program writes the data to a text file.

    PARAMETERS
    ----------
    sp2_inputs : dict
        of the inputs to be used in the C code template

    RETURNS
    -------
    None
        since for now we pass and modify the model object

    """

    #
    #> create C code
    #
    with open(C_TEMPLATE_PATH, 'r') as f:
        C_template = f.read()

    C_code = C_template.format(**sp2_inputs)
    # print(C_code)

    os.chdir(C_CODE_PATH)

    cfname = 'tmp.c'
    with open(cfname, 'w') as f:
        f.writelines(C_code)

    #
    #> compile C code
    #
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
    cmd = f"gcc {cfname:s} spectrl2.c solpos.c -o {xfname:s} -lm"
    # ^ I guess I was using a quite old gcc before (linking math was not required)
    # https://stackoverflow.com/a/5005419

    subprocess.run(cmd.split())
    os.chdir(CWD)

    #
    #> run the C code
    #> catching output in the file defined by raw_ofname
    #

    raw_ofname = "tmp.csv"

    # TODO: use tempfile from std lib to generate unique files to avoid possible conflict

    cmd2 = f"{C_CODE_PATH:s}/{xfname:s}"
    with open(raw_ofname, 'w') as f:
        subprocess.run(cmd2, stdout=f)

    #
    #> load raw data from the SPCTRAL2 run
    #> 
    #> first row in raw file is NREL SPA computed solar zenith angle
    #> the rest of the rows are three cols: wavelength, direct, diffuse
    #
    with open(raw_ofname, 'r') as f:
        sza = float(f.readline().rstrip())
    wl, Idr, Idf = np.loadtxt(raw_ofname, delimiter=',', skiprows=1, unpack=True)

    return {
        "sza": sza, "wl": wl, "Idr": Idr, "Idf": Idf
    }



class model():
    """Wrapper for the SPCTRAL2 model (C version).
    
    Note that there are no runtime options in the original model,
    so we compile every time (not very efficient...)
    but luckily it takes very little time to run. 

    Note: the init arguments are set as keyword-only and so can be passed in any order, 
    and any selection of them can be used. 

    Initialization parameters
    -------------------------
    lat, lon
        location (deg.)
        lon as within [-180, 180]
    year, month, day, hour, minute, second
        local time! (integers only)
    utcoffset
        UTC offset (negative for west)
    temp
        ambient temperature (deg. C)
    press
        pressure (surface pressure at the location) (mb/hPa)
    tau500
        aerosol optical depth at 0.5 μm (500 nm), base e
        set to 0.0 to ignore (-1.0 to silently fail with warning in output file...)
        default value (from original C example run script) = 0.2, from Excel version = 0.27 (midrange)
        typical range: [0.05, 0.55] for clear sky
    watvap
        column precipitable water vapor (cm)
        set to 0.0 to ignore (-1.0 to silently fail with warning in output file...)
        default value (from original C example run script) = 1.36, from Excel version = 1.42
        typical range: [0.3, 6]
    ozone
        total column ozone (cm)
        set to -1.0 for an internal (SPCTRAL2) estimation based on lat/lon and time of day/year (Heuklon 1978)
        default value (from original C example run script) = -1.0

    casename
        name of the case, used to create a directory for saving the outputs
        default = 'test'
    ID
        identifier used for the saved output filename
        default = '001'

    """

    def __init__(self, *, 
        lat: float = 40.0, 
        lon: float = -105.0, 
        year: int = 2014, 
        month: int = 7, 
        day: int = 1,
        hour: int = 12, 
        minute: int = 0, 
        second: int = 0,
        utcoffset: int = -5.,
        #
        temp: float = 27.0, 
        press: float = 1010.0,
        #
        tau500: float = 0.27, 
        watvap: float = 1.42, 
        ozone: float = -1.0, 
        #
        casename: str = 'test', 
        ID: str = '001',
        #
        **kwargs,
        ):

        val_wri = validate_within_range_inclusive  # alias

        # validate lat/lon
        val_wri(lat, (-90, 90), name="lat")
        val_wri(lon, (-180, 180), name="lon")

        # validate date/time inputs by trying to make a datetime
        self.dt = datetime.datetime(year, month, day, hour, minute, second)
        # TODO: could incorporate timezone info (with utcoffset)?

        # utcoffset
        val_wri(utcoffset, (-12, 12), name="utcoffset")
        # note: technically can be +14 max (https://en.wikipedia.org/wiki/List_of_UTC_time_offsets)

        # temperature (recall deg. C, supposed to be near-sfc)
        val_wri(temp, (-273.15 + 1e-6, None), (-70.0, 100.0), name="temp")

        # pressure (recall mb, supposed to be near-sfc pressure)
        val_wri(press, (1e-6, None), (500, 1200), name="press")

        # tau500 (supposed to be a clear sky model so should be low)
        val_wri(tau500, (0, None), (0.01, 1.0), name="tau500")

        # watvap
        val_wri(watvap, (0, None), (0.1, 10.0), name="watvap")
        
        # ozone is special since -1 is an allowed value as a flag
        # (that we use by default)
        # so neglecting validating for now (until update the validator to allow this)

        # since we are using the same variable names, 
        # we can grab the ones we need this way 
        current_locals = {k:v for k, v in locals().items()}
        self.sp2_inputs = {
            name: current_locals[name]
            for name in SPCTRAL2_MODEL_REQUIRED_INPUT_PARAMETER_NAMES
        }
        # ^ could be one expression, iterating over locals with an if
        # TODO: separate validation from __init__ so safely change values and run again?

        # in prep for running
        self._has_been_run = False
        self._corrected = False
        self._out_raw = False
        self.out = False        


    def run(self, *, sza_check=None):
        """ 
        Run the SPCTRAL2 C program using the selected method (C code or Cython). 

        Parameters
        ----------
        sza_check : float, optional
            SZA (solar zenith angle) value calculated by another program
            to check that NREL SOLPOS is working satisfactorily
        
        """
        
        # for now, running again not allowed
        if self._has_been_run:
            raise Exception("Create a new model object to run again.")

        # run using selected runner
        out = _recompile_run(self.sp2_inputs)  # only this one is set up

        self._has_been_run = True

        #> check sza (if check value provided)
        #  first row in raw file is solar zenith angle
        if sza_check is not None:
            assert isinstance(sza_check, (int, float))

            sza_sp2 = out["sza"]

            if abs(sza_sp2 - sza_check) > 0.1:  # should be very close
                print('SZA issue:')
                print(f'  `sp2`  : {sza_sp2:.3f}')
                print(f'  `check`: {sza_check:.3f}')
            else:
                print('SZA matches check.')

        #
        #> create output xarray with the calculated values
        #  (could be separate fn) 
        
        # collect from the SPCTRAL2 output
        wl = out["wl"]
        Idr = out["Idr"]
        Idf = out["Idf"]
        sza = out["sza"]

        time = self.dt  # datetime
        sdt = self.dt.strftime('%Y-%m-%d %H:%M:%S')

        # dwl (wavelength band widths) should be defined at wlc (wavelength band centers)
        # but this depends on integration/interpolation method used,
        # so for now we don't do it here

        # we want it to be easily concat-able (in case running for multiple times)
        # so using a time coordinate variable (with one value)
        # we have to modify the arrays to support this (adding dim)
        # -> may undo this at some point
        dset = xr.Dataset(
            coords={'wl': ('wl', wl, {'units': 'μm', 'long_name': 'wavelength'}), 
                    'time': np.atleast_1d(time)}, 
            data_vars={
                'Idr': (
                    ('time', 'wl'), Idr[np.newaxis, :], 
                    {'units': f'W m^-2 μm^-1', 'long_name': 'direct beam solar irradiance'}
                ),
                'Idf': (
                    ('time', 'wl'), Idf[np.newaxis, :], 
                    {'units': f'W m^-2 μm^-1', 'long_name': 'diffuse solar irradiance'}
                ),
                'sdt': ('time', np.atleast_1d(sdt), {'long_name': 'datetime string'}),
                'sza': ('time', np.atleast_1d(sza), {'units': 'deg.', 'long_name': 'solar zenith angle'}),
            },
            attrs={
                "sp2_inputs": self.sp2_inputs
            },
        )

        #> return and store (probably only need to do one)
        self._out_raw = dset
        self.out = self._out_raw  # until correction happens?

        # to enable chains like `model().run().out` and `model().run().plot()`
        return self


    def plot(self):
        """Plot the output spectra."""
        
        if not self._has_been_run:
            raise Exception("Gotta run the model before plotting the results.")

        plot_spectra(self.out, save=False)



def plot_spectra(dset, *,
    title='', title_left='', title_right='',
    save=False):
    """Plot spectrum from an output CSV file, given the filename.

    Parameters
    ----------
    dset : xr.Dataset
        the model output, stored in attribute `model.out`
    title, title_left, title_right : str
        additional info to be added to plot
    save : bool or str, optional
        can use to supply the figure save path (without extension)

    Returns
    -------
    pyplot figure handle
    """
    import matplotlib.pyplot as plt

    wl = dset.wl
    Idr = dset.Idr.squeeze()
    Idf = dset.Idf.squeeze()

    f1, a = plt.subplots(figsize=(8, 4))

    a.plot(wl, Idr, label='direct')
    a.plot(wl, Idf, label='diffuse')
    a.plot(wl, Idr+Idf, label='total')
    
    # note: SPCTRAL2 limits are [0.3, 4.0]
    a.set_xlim((0.3, 2.6))  

    a.set_title(title)
    a.set_title(title_left, loc='left')
    a.set_title(title_right, loc='right')
    a.set_xlabel('wavelength (μm)')
    a.set_ylabel('spectral irradiance (W m$^{-2}$ μm$^{-1}$)')
    
    a.legend()
    a.grid(True)
    f1.tight_layout()

    if save:
        raise NotImplementedError

        # if isinstance(save, str):
        #     try:
        #         p = Path(save)
        #     except:

        # outdir = './out/{casename:s}/img'.format(casename=casename)
        # f1.savefig('{:s}/{:s}.png'.format(outdir, ID), dpi=150, 
        #     transparent=False, bbox_inches='tight', pad_inches=0.05)

    return f1



def correct(dset, 
    measured_units='E', measured_val=None, measured_bnds=(0.4, 0.7),
    total_solar=None,
    plot=False, save_plot=False,
    cc_correct=False, cc_p1=1.0, cc_p2=0.3,
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
    dset : xr.Dataset

    measured_units : str
        'E' (W/m^2) or 'photons' (photon flux density)
    measured_val : float
        A measured broadband irradiance value (μmol photons / m^2 / s)
    measured_bnds : tuple
        Wavelength band bounds for the measured irradiance value (μm) 
    total_solar : float
        A measurement of total solar irradiance (direct+diffuse; W/m^2)
        We can use this to correct the regions outside the measured band (e.g., PAR)
        Only use this in addition to a measured band value!
    cc_correct : bool
        Whether to apply the (experimental) cloud corrections
    cc_p1 : float
        Bird et al. (1987) cloud-correction parameter 1
        for wavelengths in range <= 0.55 μm
    cc_p1 : float
        Bird et al. (1987) cloud-correction parameter 2
        for wavelengths in range 0.5--0.926 μm

    """
    # TODO: move cloud correction to separate fn (can still pass kwargs to it from this one)
    # TODO: cloud correction fn should optionally run with the watvap enhancement previously in __init__:
    # if self.cc_correct:
    #     self.watvap = watvap * 1.2 * 1.2 * self.cc_p1
    #     #^ climatologically higher watvap for cloudy conditions + 
    #     #  increase in water vapor path due to multiple scattering within cloud

    #files = glob('{:s}/raw/*.csv'.format(self.outdir))

    #> load raw data from the SPCTRAL2 run 
    #  first row in raw file is NREL SPA computed solar zenith angle
#        with open(self.raw_ofname, 'r') as f:
#            sza = f.readline().rstrip()
    # wl, Idr, Idf = np.loadtxt(self.raw_ofname, delimiter=',', skiprows=1, unpack=True)

    wl = dset.wl
    Idr = dset.Idr.squeeze()
    Idf = dset.Idf.squeeze()

    #> NaN -> 0
    # TODO: check if this is really necessary
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
    if cc_correct:
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
        Idf += Idr * (1-cc_p2)
        Idr -= Idr * (1-cc_p2)
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
        # convert SPCTRAL2 sub-band irradiances in W/m^2/μm to broadband in photon flux units (photons/m^2/s)
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

    # #> save corrected spectra
    # fnew = self.corrected_ofname
    # a = np.vstack((wl, dwl, Idr2, Idf2)).T
    # #print a.shape
    # np.savetxt(fnew, a, delimiter=',',
    #     fmt=['%.3e', '%.3e', '%.6e', '%.6e'])

    # #> plot if desired
    # if plot == True:
    #     plot_spectrum(casename=self.casename, ID=self.ID, output_type='corrected',
    #         title='', title_left='', title_right='',
    #         save=save_plot)

