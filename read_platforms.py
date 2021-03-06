import numpy as np
import pandas as pd
import datetime as dt
from datetime import timedelta
import metpy.calc as mpcalc
from metpy.units import units
import xarray as xr
import glob
import warnings
warnings.filterwarnings(action='ignore')

def read_unlmm(file,tstart=None, tend=None):
    mmfile= glob.glob(file)
    print(np.shape(mmfile))
    print("mmfile = ", mmfile)

    ds = xr.open_dataset(mmfile[0])

    #convert form epoch time to utc datetime object
    timearray = np.array([dt.datetime.utcfromtimestamp(t/1e9) for t in ds.time.values])
    U,V = mpcalc.wind_components(ds.wind_speed, ds.wind_dir)

    lats = np.array([40])*units.degrees
    #create new xarray dataset to make plotting easier cause original is trash
    dims = ['datetime']
    coords = {
        'datetime': timearray
    }
    data_vars = {
        'lat': (dims, ds.lat.values, {'units':str(lats.units)}),
        'lon': (dims, ds.lon.values, {'units':str(lats.units)}),
        'Z_ASL': (dims, ds.alt.values, {'units':str(ds.alt.units)}),
        'Z_AGL': (dims, np.zeros_like(ds.alt.values), {'units':str(ds.alt.units)}),
        'Temperature': (dims, ds.fast_temp.values, {'units':str(ds.fast_temp.units)}),
        'Dewpoint': (dims, ds.dewpoint.values, {'units':str(ds.dewpoint.units)}),
        'RH': (dims, ds.calc_corr_RH.values, {'units':str(ds.calc_corr_RH.units)}),
        'Pressure': (dims, ds.pressure.values, {'units':str(ds.pressure.units)}),
        'U': (dims, U.m, {'units':str(U.units)}),
        'V': (dims, V.m, {'units':str(V.units)}),
        'Theta': (dims, ds.theta.values, {'units':str(ds.theta.units)}),
        'Thetav': (dims, ds.theta_v.values, {'units':str(ds.theta_v.units)}),
        'Thetae': (dims, ds.theta_e.values, {'units':str(ds.theta_e.units)})
    }
    subds = xr.Dataset(data_vars, coords)

    #convert to pandas
    pd_unl=subds.to_dataframe()
    pd_unl.reset_index(inplace=True)

    # Subset the dataset to the desired time
    if tstart is None:
        hstart = pd_unl['datetime'].iloc[0]
    else:
        hstart = tstart

    if tend is None:
        hend= pd_unl['datetime'].iloc[-1]
    else:
        hend   = tend

    # Save only desired iop data
    data_unl = pd_unl.loc[(pd_unl['datetime']>=hstart)&(pd_unl['datetime']<=hend)]

    return data_unl

# *******************

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
    data_nssl = data.loc[(data['time']>=hstart_dec)&(data['time']<=hend)]

    # Convert time into timedeltas
    date = dt.datetime.strptime('2019-'+file[-15:-13]+'-'+file[-13:-11], '%Y-%m-%d')
    time_deltas = []
    for i in np.arange(len(data_nssl)):
        j = data_nssl['time'].iloc[i]
        time_deltas = np.append(time_deltas,date+dt.timedelta(hours=int(j),minutes=int((j*60) % 60),seconds=int((j*3600) % 60)))
    data_nssl.loc[:,'datetime']=time_deltas  # This form is faster and prevents warnings

    # Caclulate desired variables
    p,t   = data_nssl['p'].values*units.hectopascal, data_nssl['tfast'].values*units.degC
    theta = mpcalc.potential_temperature(p,t)
    data_nssl.loc[:,'Theta'] = theta.magnitude

    r_h     = data_nssl['rh'].values/100
    mixing = mpcalc.mixing_ratio_from_relative_humidity(r_h,t,p)
    thetav = mpcalc.virtual_potential_temperature(p,t,mixing)
    data_nssl.loc[:,'Thetav'] = thetav.magnitude

    td     = mpcalc.dewpoint_rh(temperature=t, rh=r_h)
    thetae = mpcalc.equivalent_potential_temperature(p,t,td)
    data_nssl.loc[:,'Thetae'] = thetae.magnitude

    Spd =data_nssl['spd'].values*units('m/s')
    dire =data_nssl['dir'].values*units('degrees')
    u,v =mpcalc.wind_components(Spd,dire)
    #data_nssl['U']=u.to('knot')
    data_nssl.loc[:,'U'] = u.to('knot')

    #data_nssl['V']=v.to('knot')
    data_nssl.loc[:,'V'] = v.to('knot')

    q_list=['qc1','qc2','qc3','qc4']
    data_nssl.loc[:,'qc_flag']=data_nssl[q_list].sum(axis=1)

    return data_nssl

# *******************

def maskdata(p_var, platform_file, mask=True):
    ''' Read in a dataset and return the masked version of it
            Masking is based of QC flags etc and can be diff for each platform
            If mask= False then the defn will simply return the dataset unmodified
    '''
    platform_unmasked= platform_file[p_var].values

    platform_data = platform_unmasked

    if mask== False:
        platform_data= platform_unmasked
    elif mask== True:
        name = platform_file.get('name', 'fixmenow3')
        platform_name= str(name)
        print("platform_name = ", platform_name)

        if (platform_name in ['FFld_df','WinS_df','LIDR_df','Prb1_df','Prb2_df']):
            platform_data= np.ma.masked_where(platform_file['qc_flag'].values>0, platform_unmasked)

        # This was needed to work with the nextrad code in it's present form
        elif (platform_name in ['FFld','WinS','LIDR','Prb1','Prb2']):
            platform_data= np.ma.masked_where(platform_file['qc_flag'].values>0, platform_unmasked)

        elif (platform_name in ['CoMeT1_df','CoMeT2_df','CoMeT3_df']):
            #for now
            platform_data= platform_unmasked

        # This was needed to work with the nextrad code in it's present form
        elif (platform_name in ['CoMeT1','CoMeT2','CoMeT3']):
            #for now
            platform_data= platform_unmasked

        elif (platform_name in ['Insert UAS filenames here']):
            print("will be filled in")
        else:
            print("What platform are you trying to use?")

    return platform_data
