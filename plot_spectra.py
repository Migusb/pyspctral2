
"""

plot each of the corrected spectra to make sure they make sense
for CR 201407 


"""

from glob import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

plt.close('all')

micron = u'\u03BCm'  

lw = 1.5

ylim_max = 1500

def plot1(fpath, title='', outdir='./'):

    fbnwext = os.path.basename(fpath)
    fbn, fext = os.path.splitext(fbnwext)
    ID = fbn[-3:]

    wl, dwl, direct, diffuse = np.loadtxt(fpath, delimiter=',', skiprows=0, unpack=True)



    plt.figure(figsize=(6,4.0))

    plt.plot(wl, direct+diffuse, lw=lw, label='total')
    plt.plot(wl, direct, '--', lw=lw, label='direct')
    plt.plot(wl, diffuse, ':', lw=lw, label='diffuse')

    #plt.xlabel(r'wavelength, \si{\micro\metre}')
    plt.xlabel(r'wavelength ('+micron+')')
    #plt.ylabel(r'irradiance, \si{W.m^{-2}.\micro\meter^{-1}}')
    plt.ylabel(r'irradiance (W m$^{-2}$ '+micron+')')

    plt.title(title)

    plt.xlim(0.30, 2.6)  # same bounds as idealized leaf rad props data, now extended to 0.30
    #plt.ylim(0, plt.gca().get_ylim()[1])
    plt.ylim(0, ylim_max)

    plt.legend()

    plt.tight_layout()

    plt.savefig('{:s}/{:s}.png'.format(outdir, ID), dpi=150, 
                transparent=False, bbox_inches='tight', pad_inches=0.05)

    plt.close('all')



env_data_f = '/storage/home/zlm1/work/access/zm/data/CH14_186-188.dat'

sdate, stime = np.genfromtxt(env_data_f, dtype=str,
    skip_header=15, skip_footer=0,
    usecols=(0, 1), unpack=True,
    )

outdir = './CR_201407_3day/img'

files = glob('./CR_201407_3day/corrected/*.csv')

for i, f in tqdm(enumerate(sorted(files))):  # enumerate messes up tqdm bar

    s = sdate[i] + ' ' + stime[i]
    plot1(f, title=s, outdir=outdir)










