import numpy as np
import pandas as pd
import datetime as dt
from datetime import timedelta
import metpy.calc as mpcalc
from metpy.units import units
import warnings
warnings.filterwarnings(action='ignore')

def read_nsslmm(file,tstart=None,tend=None):
    ''' Reads data files provided by NSSL on TORUS19 EOL site
        Important that datafiles and filenames follow TORUS19 readme
        INPUT: filename string (following NSSL naming convention)
               tstart,tend: datetime objects
        OUTPUT: dataframe containing all data collected during desired times
    '''
    column_names = ['id','time','lat','lon','alt',
                    'tfast','tslow','rh','p','dir',
                    'spd','qc1','qc2','qc3','qc4']
    # Read NSSL file using column names from readme
    data = pd.read_csv(file,header=0,delim_whitespace=True,names=column_names)
    data = data.drop_duplicates()

    # Find timedelta of hours since start of iop
    # (IOP date taken from filename!)
    tiop = dt.datetime(2019, np.int(file[-15:-13]),np.int(file[-13:-11]),0,0,0)
    if tstart is None:
        hstart = tiop
        #convert to decimal hours HH.HHH
        hstart_dec=hstart.hour+(hstart.minute/60)+(hstart.second/3600)
    else:
        hstart = (tstart-tiop).seconds/3600
        hstart_dec=hstart

    if tend is None:
        hend = data['time'].iloc[-1]
    else:
        hend   = (tend-tiop)
        if hend >= dt.timedelta(days=1):
            hend = (tend-tiop).seconds/3600 + 24.
        else:
            hend = (tend-tiop).seconds/3600
    # Save only desired iop data
    data_iop = data.loc[(data['time']>=hstart_dec)&(data['time']<=hend)]

    # Convert time into timedeltas
    date = dt.datetime.strptime('2019-'+file[-15:-13]+'-'+file[-13:-11], '%Y-%m-%d')
    time_deltas = []
    for i in np.arange(len(data_iop)):
        j = data_iop['time'].iloc[i]
        time_deltas = np.append(time_deltas,date+dt.timedelta(hours=int(j),minutes=int((j*60) % 60),seconds=int((j*3600) % 60)))
    data_iop['datetime']=time_deltas

    # Caclulate desired variables
    p,t   = data_iop['p'].values*units.hectopascal, data_iop['tfast'].values*units.degC
    theta = mpcalc.potential_temperature(p,t)
    data_iop['Theta'] = theta.magnitude

    r_h     = data_iop['rh'].values/100
    mixing = mpcalc.mixing_ratio_from_relative_humidity(r_h,t,p)
    thetav = mpcalc.virtual_potential_temperature(p,t,mixing)
    data_iop['Thetav'] = thetav.magnitude

    td     = mpcalc.dewpoint_rh(temperature=t, rh=r_h)
    thetae = mpcalc.equivalent_potential_temperature(p,t,td)
    data_iop['Thetae'] = thetae.magnitude

    Spd =data_iop['spd'].values*units('m/s')
    dire =data_iop['dir'].values*units('degrees')
    u,v =mpcalc.wind_components(Spd,dire)
    data_iop['U']=u.to('knot')
    data_iop['V']=v.to('knot')

    q_list=['qc1','qc2','qc3','qc4']
    data_iop['qc_flag']=data_iop[q_list].sum(axis=1)

    return data_iop

