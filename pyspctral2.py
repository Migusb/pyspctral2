
"""

Python module to run the SPCTRAL2 C code and catch output

"""

import datetime as dt
from glob import glob
import os
from pathlib2 import Path
from subprocess import call

import numpy as np
import xarray as xr


C_code_path = './C_code'  # loc of SPCTRAL2 and solpos C code
cwd = os.getcwd()

micron = u'\u03BCm'

with open('./C_code/run_spctral2_template.c', 'r') as f:
    C_template = f.read()


class model():
    """ Wrapper for the SPCTRAL2 model (C version)
    
    Note that there are no runtime options in the original model,
    so we compile every time (not very efficient...)

    """

    def __init__(self, 
        lat=45.0, lon=-80.0, 
        year=2014, month=7, day=1,
        hour=12, minute=0, second=0,
        utcoffset=-5.,
        temp=20, press=1000,
        casename='test', ID='001'
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
 
        self.temp = temp  # near-sfc air temp (deg. C)
        self.press = press  # near-sfc air pressure (mb)

        # ozone, water, etc

        #> get ready for output

        self.casename = casename  # group of runs
        self.ID = ID  # identifier for one run of the group
        dirbases = ['img', 'raw', 'corrected']
        dirroot = './out/{:s}'.format(casename)

        self.outdir = dirroot

        for dirbase in dirbases:
            #os.makedirs(dirroot+dirbase, exist_ok=True)  # exist_ok only in py3 version
            ps = '{:s}/{:s}'.format(dirroot, dirbase)
            Path(ps).mkdir(parents=True, exist_ok=True)

        self.raw_ofname = '{:s}/raw/{:s}.csv'.format(self.outdir, self.ID)
        self.corrected_ofname = '{:s}/corrected/{:s}.csv'.format(self.outdir, self.ID)
        


    def run(self):
        """ """
        
        #> create C code
        C_code = C_template.format(lon=self.lon, lat=self.lat,
                year=self.year, month=self.month, day=self.day,
                hour=self.hour, minute=self.minute, second=self.second,
                utcoffset=self.utcoffset,
                temp=self.temp, press=self.press)
        #print C_code
        
        os.chdir(C_code_path)

        cfname = 'tmp.c'
        with open(cfname, 'w') as f:
            f.writelines(C_code)

        #> compile C code
        xfname = 'tmp.out'
        cmd = 'icc {cfname:s} spectrl2.c spectrl2.h solpos.c solpos00.h -o {xfname:s}'.format(cfname=cfname, xfname=xfname) 
        call(cmd.split())
        os.chdir(cwd)

        #> run the C code
        #  catching output in ofname

        #cmd2 = './{xfname:s} > {ofname:s}'.format(xfname=xfname, ofname=ofname)
        cmd2 = '{c:s}/{xfname:s}'.format(xfname=xfname, c=C_code_path)
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
            Assumes units umol photons / m^2 / s for the measured value


        """

        #files = glob('{:s}/raw/*.csv'.format(self.outdir))


        # first row in raw file is solar zenith angle
        with open(self.raw_ofname, 'r') as f:
            sza = f.readline().rstrip()
        wl, Idr, Idf = np.loadtxt(self.raw_ofname, delimiter=',', skiprows=1, unpack=True)

        #print sza

        Idr[np.isnan(Idr)] = 0.0
        Idf[np.isnan(Idf)] = 0.0

        # correct for unrealistic and unphysical values
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



        dwl = np.diff(wl)
        dwl = np.append(dwl, dwl[-1])  # or np.nan

        assert( len(measured_bnds) == 2 and measured_bnds[0] < measured_bnds[1] )
        wl_meas = (wl >= measured_bnds[0]) & (wl <= measured_bnds[1])

        I_meas_sp2_all = (Idr+Idf)[wl_meas]  # all bands within the measured broad band
        dwl_meas = dwl[wl_meas]

        I_meas_sp2 = ( I_meas_sp2_all*dwl_meas / (6.626e-34*2.998e8/(wl[wl_meas]*1e-6)) ).sum() / 6.022e23 * 1e6

        cf = measured_val / I_meas_sp2  # correction factor using measured broadband irradiance value

        Idr2 = Idr * cf
        Idf2 = Idf * cf 

        #> recheck the broadband calculation with sp2 data to make sure correction worked

        I_meas_sp2_all_2 = (Idr2+Idf2)[wl_meas]

        I_meas_sp2_2 = ( I_meas_sp2_all_2*dwl_meas / (6.626e-34*2.998e8/(wl[wl_meas]*1e-6)) ).sum() / 6.022e23 * 1e6

        assert( np.all(np.isclose(I_meas_sp2_2, measured_val)) )


        #> save corrected spectra
        
        fnew = self.corrected_ofname

        a = np.vstack((wl, dwl, Idr2, Idf2)).T
        #print a.shape
        np.savetxt(fnew, a, delimiter=',',
            fmt=['%.3e', '%.3e', '%.6e', '%.6e'])


        #> plot if desired
        #if plot == True:
        
            
def combine(casename='blah', time=[], info=''):
    """
    For given casename, combine all of the individual output files and data
    into one netCDF file

    time: assumed list of datetime objs

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
            szai = f.readline().rstrip()
        sza.append(szai)

    #> build Idr and Idf arrays
    Idr = np.zeros((nt, nwl))
    Idf = np.zeros_like(Idr)
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
            },
        attrs={'casename': casename, 'info': info},
        )

    #> save dataset
    ofname = './out/{:s}.nc'.format(casename)
    dset.to_netcdf(ofname)



