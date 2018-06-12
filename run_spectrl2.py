
"""

run the C version of SPCTRAL2 for each of the env data time steps 
for Chestnut Ridge 5-7 July sim


"""

import datetime as dt
#from glob import glob
#import os
#from subprocess import call

import numpy as np
from tqdm import tqdm

import pyspctral2 as psp2
reload(psp2)  # < so changes take effect when the module is edited..

#> location of Chesnut Ridge canopy site
lon = -84.2875
lat = 35.9583

#> construct datetimes manually
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

#def run_spctral2_all():


casename = 'CR_201407_3day'

dts = []

for i, line in enumerate(tqdm(env_data[:])):
    
    sdate, stime, Tair_C, p_mb, PPFD = line  # unpack

    ID = '{:03d}'.format(i)

    sdt = sdate + ' ' + stime
    dtline = dt.datetime.strptime(sdt, '%Y-%m-%d %H:%M:%S')
    #print dtline
    dts.append(dtline)

    year   = dtline.year
    month  = dtline.month
    day    = dtline.day
    hour   = dtline.hour
    minute = dtline.minute
    second = dtline.second

    #print year, month, day, hour, minute, second
 

#    model = psp2.model(lon=lon, lat=lat,
#            year=year, month=month, day=day,
#            hour=hour, minute=minute, second=second,
#            utcoffset=-5.,
#            temp=Tair_C, press=p_mb,
#            casename=casename, ID=ID)
#
    #model.run()

    #model.correct(measured_val=PPFD, measured_bnds=(0.3, 0.7))


#> combine

psp2.combine(casename=casename, time=dts, info='info.')


#if __name__ == '__main__':

    #run_spctral2_all()


