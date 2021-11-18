#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import needed modules
######################
import matplotlib
#  matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.patheffects as PathEffects
from matplotlib import ticker
from matplotlib.ticker import (LinearLocator, FixedLocator, MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.cm as cmx
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler
import cartopy
import time
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import ShapelyFeature, NaturalEarthFeature
from cartopy.io.shapereader import Reader
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
import metpy
from metpy.plots import StationPlot, ctables
from metpy.units import units
import metpy.calc as mpcalc
#  import modin.pandas as pd
import pandas as pd
from pandas.plotting import register_matplotlib_converters
import numpy as np
import xarray as xr
#from netCDF4 import num2date
import os
import os.path
import os
import glob
from pathlib import Path
from joblib import Memory, Parallel, delayed
import scipy
from scipy import ndimage, interpolate
import cmocean
import gc
import nexradaws
import pyart
import sys
import cftime # for num2pydate
import pprint
import math
import copy
import fnmatch
#  import pytda
#  import singledop
import fire
#  from PyMeso.pymeso import llsd
from astropy.convolution import convolve

pp = pprint.PrettyPrinter(indent=4)
import numpy as np
import matplotlib.pyplot as plt
import pyart
from netCDF4 import Dataset
from metpy.plots import ctables
import cmocean
import os
import matplotlib.patheffects as PathEffects
import pandas as pd

import warnings
import llsd
warnings.filterwarnings("ignore")

from plotwdssiiQuick import readwdssii

#
def radar_fields_prep(config, rfile, radar_type, sweep_id, moment='From_config'):
    ###########
    #  print(rfile.info(level='full'))
    #  print(rfile.info(level='standard'))
    #  print(rfile.info(level='compact'))

    if radar_type == 'KA':
        vel_name, specw_name = 'corrected_velocity','spectrum_width'
        ncp_name, ncp_thresh = 'normalized_coherent_power', .4
        refl_name, refl_thresh = 'reflectivity', -5

    if radar_type == 'NOXP':
        vel_name , refl_name, specw_name = 'corrected_velocity', 'DBZ','WIDTH'
        ncp_name, ncp_thresh = 'SQI', .6
        refl_name, refl_thresh = 'DBZ', 15
    if radar_type == 'WSR':
        if config.radar_controls['Use_downloaded_files']== True:
            vel_name , refl_name, specw_name = 'corrected_velocity', 'reflectivity','spectrum_width'
        else:
            vel_name , refl_name, specw_name = 'velocity', 'reflectivity','spectrum_width'
        refl_thresh =15
    #creating the mask for attenuation
        #  range_mask = np.zeros(np.shape(reflectivity))
        #  for i in range(0, len(range_mask[:,0])): range_mask[i,:] = rfile.range['data'] > (rfile.range['data'][-1]-1000.)
        #  range_mask = range_mask.astype(bool)
        #  total_mask = [any(t) for t in zip(range_mask.flatten(), normal_mask.flatten())]
        #  gatefilter = pyart.filters.GateFilter(rfile)
        #  gatefilter.exclude_below('DBZ', 10)
    #  gatefilter= pyart.filters.moment_based_gate_filter(rfile, ncp_field= ncp_name, rhv_field=rhv_name,
                                                      #  refl_field=refl_name, min_ncp= ncp_thresh)

    def mask(rfile, gatefilter, moment, new_mom_name, masking=True):
        if masking==True:
            #  gatefilter=pyart.correct.despeckle_field(rfile, moment, gatefilter=gatefilter, size= 15)
            mom_masked= np.ma.masked_where(gatefilter.gate_included == False, moment)
        else: mom_masked=moment
        rfile.add_field(new_mom_name, {'data': mom_masked}, replace_existing=True)
        return rfile

    gatefilter = pyart.filters.GateFilter(rfile)
    texture_feild= pyart.retrieve.calculate_velocity_texture(rfile, vel_field=vel_name, wind_size=3)
    rfile=mask(rfile, gatefilter, texture_feild['data'], 'vel_texture', masking=False)
    gatefilter.exclude_above('vel_texture', 3)
    gatefilter.include_above(refl_name, refl_thresh)

    vel_field = rfile.fields[vel_name]['data']
    smooth_data = scipy.ndimage.filters.median_filter(vel_field, 3)
    #  smooth_data = scipy.ndimage.filters.median_filter(vel_field, 3)
    #  smooth_data_ma = np.ma.masked_where(np.ma.getmask(vel_field), smooth_data)
    rfile=mask(rfile, gatefilter, smooth_data, 'vel_smooth', masking =True)
#
    if moment == 'From_config': mom = config.r_mom
    else: mom = [moment]
    for m in mom:
        print(m)

        if m in ['vel', 'refl', 'specw']:
            if m == 'vel':
                orig_field_name, new_feild_name, gfilter = vel_name, 'vel_fix', gatefilter
            if m == 'refl':
                if radar_type != 'WSR':
                    gatefilter_refl= pyart.filters.moment_based_gate_filter(rfile, ncp_field= ncp_name, min_ncp= ncp_thresh)
                else:
                    gatefilter_refl=gatefilter
                orig_field_name, new_feild_name, gfilter = refl_name, 'refl_fix', gatefilter_refl
            if m == 'specw': orig_field_name, new_feild_name, gfilter = specw_name, 'specw_fix', gatefilter
            orig_field = rfile.fields[orig_field_name]['data']
            rfile=mask(rfile, gfilter, orig_field, new_feild_name, masking=True)

        if m in ['vel_despeck', 'refl_despeck']:
            if m == 'refl_despeck': orig_field_name, new_feild_name= 'refl_fix', 'refl_despeck'
            elif m == 'vel_despeck': orig_field_name, new_feild_name= 'vel_fix', 'vel_despeck'
            speckeled_gatefilter = pyart.correct.despeckle_field(rfile, field=orig_field_name, threshold=-100,
                                                                     gatefilter=gatefilter, size=30)
            rfile=mask(rfile, speckeled_gatefilter, rfile.fields[orig_field_name]['data'], 'refl_despeck')

        # * * *
        if m =='vel_test':
            #  data_masktest= np.ma.masked_array(rfile.fields[vel_name]['data'], gatefilter.gate_excluded)
            #  rfile.add_field('vel_gfilt', {'data': data_masktest}, replace_existing=True)
            #  texture_feild= pyart.retrieve.calculate_velocity_texture(rfile, vel_field='vel_gfilt', wind_size=3)
            #  texture_feild= pyart.retrieve.calculate_velocity_texture(rfile, vel_field='vel_savgol', wind_size=3)
            texture_feild= pyart.retrieve.calculate_velocity_texture(rfile, vel_field=vel_name, wind_size=3)
            rfile=mask(rfile, gatefilter, texture_feild['data'], m, masking=True)

        # * * *
        if m == 'vel_texture_dea':
            texture_feild= pyart.retrieve.calculate_velocity_texture(rfile, vel_field=vel_name, wind_size=3)
            rfile=mask(rfile, gatefilter, texture_feild['data'], m, masking=False)
        if m == 'vel_texture_smoothed':
            #  vel_field = rfile.fields['vel_fix']['data']
            #  vel_savg = scipy.signal.savgol_filter(vel_field, window_length=13, polyorder=2, axis=1)
            #  rfile=mask(rfile, gatefilter, vel_savg, 'vel_savgol')
            texture_feild= pyart.retrieve.calculate_velocity_texture(rfile, vel_field='vel_smooth', wind_size=3)
            smooth_data = scipy.ndimage.filters.median_filter(texture_feild['data'], 3)
            rfile=mask(rfile, gatefilter, smooth_data, m, masking=False)

        # * * *
        if m == 'vel_grad':
            vel_field = rfile.fields['corrected_velocity']['data']
            vel_smoothed=aliasfix(vel_field,13,rfile.get_nyquist_vel(0))
            vel_grad= np.gradient(vel_smoothed, 15, axis= 1)*100
            rfile=mask(rfile, gatefilter, vel_grad, 'vel_gradient')

        if m in ['vel_savgol', 'vel_savgol_axis0', 'vel_savgol_axis1']:
            if m == 'vel_savgol':
                orig_field_name = 'vel_fix'
                DER = 0
            else:
                orig_field_name = 'corrected_velocity'
                DER = 1
                if m == 'vel_savgol_axis0': A = 0
                if m == 'vel_savgol_axis1': A = 1

            vel_field = rfile.fields[orig_field_name]['data']
            PO, WL = 2, 13
            vel_savg = scipy.signal.savgol_filter(vel_field, window_length=WL, polyorder=PO, derivative=DER, axis=A)
            rfile = mask(rfile, gatefilter, vel_savg, m)

        # * * *
        if m == 'sim_vel':
            sd_test = singledop.SingleDoppler2D(L=30.0, radar=rfile, range_limits=[0, 20], sweep_number=0,
                                                name_vr='velocity', thin_factor=[4, 12])
            pyart.util.simulated_vel_from_profile(rfile, rfile.fields['vel_fix']['data'])
            rfile=mask(rfile, gatefilter, sd_test, m)

        # * * *
        if m in ['vort', 'vort_smooth']:
            if m == 'vort':
                vel_field= vel_name
            if m == 'vort_smooth':
                vel_field = rfile.fields['corrected_velocity']['data']
                smooth_data = scipy.ndimage.filters.median_filter(vel_field, 3)
                #  np.ma.masked_where(np.ma.getmask(vel_field), smooth_data)
                rfile=mask(rfile, gatefilter, smooth_data, 'vel_smooth', masking =False)
                vel_field = 'vel_smooth'

            vort_data=vort(rfile,sweep_id, vel_field)
            rfile=mask(rfile, gatefilter, vort_data, m)

        if m == 'Meso_azi':
            start_time = time.time()
            az_shear_meta = llsd.main(rfile,refl_name,vel_name)
            print("LLSD COMPUTE --- %s seconds ---" % (time.time() - start_time))
            rfile.add_field('Meso_azi', az_shear_meta, replace_existing=True)

        if m == 'az_shear':
            vel_field = rfile.fields['corrected_velocity']['data']
            smooth_data = scipy.ndimage.filters.median_filter(vel_field, 3)
            np.ma.masked_where(np.ma.getmask(vel_field), smooth_data)
            rfile=mask(rfile, gatefilter, smooth_data, 'vel_smooth', masking =True)

            az_shear_meta = llsdmain(rfile, refl_name, 'vel_smooth', sweep_id)
            rfile=mask(rfile, gatefilter, az_shear_meta['data'], m, masking=False)

        # * * *
        if m == 'diff':
            final_field= rfile.fields['vel_fix']['data']
            inital_field = rfile.fields['VEL']['data']
            diff = final_field - inital_field
            rfile=mask(rfile, gatefilter, diff, 'difference')

    return rfile

    #  gatex, gatey, gatez = radar.get_gate_x_y_z(0)
    #  rang = np.sqrt(gatex**2+gatey**2)
    #  naz,nr=np.shape(rang)
#
    #  height = gatez
#
    #  r_test = []
    #  vort = np.zeros_like(gatex,dtype=np.float64)
    #  shr = np.zeros_like(gatex,dtype=np.float64)
    #  rs = np.zeros_like(rang)

def MRMS_Reader(config, Scan_time):
    path = '/home/chris/codes/data/20190608/merged/MergedAzShear_0-2kmAGL/00.00/*'
    r_files_path= sorted(glob.glob(path))
    timeholder=[]
    for rfile in r_files_path:
        print(rfile[-13:-9])
        timeholder.append(int(rfile[-13:-9]))

    print('9999999999999')
    print(Scan_time)
    print(Scan_time.strftime('%H%M'))
    Orig_scantime=int(Scan_time.strftime('%H%M'))
    #  print(r_files_path[:, -13:-7])
    #  read_from_radar_file()
    absolute_difference_function=lambda list_value: abs(list_value-Orig_scantime)

    closest_value=min(timeholder, key=absolute_difference_function)
    #  print(str(closest_value).zfill(4))
    newtime=closest_value
    #  newtime=(str(closest_value).zfill(4))
    test=timeholder.index(newtime)
    #  print(timeholder)
    print(test)
    print(r_files_path[test])
    xi,ulat,yi,ulon,varname,refi=readwdssii(r_files_path[test])
    #  radar = pyart.io.read_cfradial(r_files_path[test], exclude_fields='time')
    #  radar = pyart.io.read(r_files_path[46], exclude_fields=['time'])

def vort(radar, swp_id, vel_name):
    '''
    Function calculates the vorticity on an PPI assuming solid body rotation. Vorticity vector is also assumed normal to plane.
    Calculation: (1/r)(1/sin(phi))*(du_r/dtheta), taken from wolfram mathworld on spherical vorticity calculations.
    Originally was negative, but due to polar coordiantes being reversed in radar world, it is actually positive.

    Parameter: radar (pyart object), swp_id (int)
    Returns: vorticity array, same shape as velocity array

    The result is multiplied by 2 due to assumed solid body rotation. If this is a bad assumption, divide the result by 2

    https://mathworld.wolfram.com/SphericalCoordinates.html
    '''
    #define the indices for the required sweep
    sweep_startidx = np.int64(radar.sweep_start_ray_index['data'][swp_id])
    sweep_endidx = np.int64(radar.sweep_end_ray_index['data'][swp_id])
    #extract data

    vel = radar.get_field(swp_id,vel_name)#in m/s
    phi = np.unique(radar.get_elevation(swp_id))[0]#in degrees
    theta = radar.get_azimuth(swp_id)#in degrees

    rangearray = np.tile(radar.range['data'],(len(theta),1))# in m
    vort_tilt=(1./(rangearray*np.sin(np.deg2rad(90-phi))))*np.gradient(vel,np.deg2rad(theta),axis=0)*2

    vort_full = np.zeros(radar.fields[vel_name]['data'].shape)
    vort_full[sweep_startidx:sweep_endidx+1] = vort_tilt
    return vort_full
'''
def ppi_vort(radar,swp_id):
    #### calculate vorticity #####
    gatex, gatey, gatez = radar.get_gate_x_y_z(0)
    rang = np.sqrt(gatex**2+gatey**2)
    naz,nr=np.shape(rang)

    height = gatez

    r_test = []
    vort = np.zeros_like(gatex,dtype=np.float64)
    shr = np.zeros_like(gatex,dtype=np.float64)
    rs = np.zeros_like(rang)

    vel = radar.fields['velocity']['data']

    for d in range(nr):#iterates thorugh all ranges
        for az in range(naz):#iterates through all azimuths
            #one with varying range (can't go past max range or beyond the radar)
               #  if abs((nr/2.)-abs(d-(nr/2.)))-1 > 0:
                   #  r2=int(abs((nr/2.)-abs(d-(nr/2.)))-1)
               #  else: r2=0
                  #sets the max size of the solid body radius to 10 bins, shouldn't need more than that, honestly could probably go smaller and save calc time
               #  if r2>10:  r=10
               #  else:  r=int(r2)
               #  rs[az,d]=2
                r = 2
                rotarr=np.zeros((nr//2),dtype=np.float64)#array of 'vorticity' from all radii
                sh = np.zeros((nr//2),dtype=np.float64)

                for i in range(r):#iterates through all possible radii
                    x1, y1= gatex[az-i,d], gatey[az-i,d]
                    x2, y2 = gatex[(az+i)%719,d], gatey[(az+i)%719,d]
                    dist = ((x2 - x1)**2 + (y2 - y1)**2)**0.5   #finds straight line dist. (m) of all possible radii

                    shear=vel[(az+i)%719,d]-vel[az-i,d]
                    rotarr[i]=(shear/dist)
                    sh[i]=(shear)
                if abs(np.nanmax(rotarr))>abs(np.nanmin(rotarr)):
                    vort[az,d]=np.nanmax(rotarr)
                    shr[az,d]=np.nanmax(sh)
                elif abs(np.nanmax(rotarr))<abs(np.nanmin(rotarr)):
                    vort[az,d]=np.nanmin(rotarr)
                    shr[az,d]=np.nanmin(sh)

    return vort

'''
def RHIvort(radar,swp_id):
    gatex, gatey, gatez = radar.get_gate_x_y_z(swp_id)
    dist = np.sqrt(gatex**2+gatey**2)
    height = gatez

    nz,nd=np.shape(dist)
    r=abs((nz/2.)-abs(np.arange(nz)-(nz/2.)))-1
    r[r<0]=0
    r=np.int_(r)#because astype didn't work... converts array to int

    vort = np.zeros_like(gatex,dtype=np.float64)
    velvort = radar.get_field(swp_id,'vel_fix')

    for d in range(nd):
        for h in range(nz):
            rotarr=np.zeros((nz//2),dtype=np.float64)
            for i in range(r[h]):
                rad=np.sqrt((dist[h+i,d]-dist[h-i,d])**2.+(height[h+i,d]-height[h-i,d])**2.)/2. #in meters
                shear=velvort[h+i,d]-velvort[h-i,d]
                rotarr[i]=np.around(shear/rad,decimals=3)

            if abs(np.nanmax(rotarr))>abs(np.nanmin(rotarr)):
                vort[h,d]=np.nanmax(rotarr)
            elif abs(np.nanmax(rotarr))<abs(np.nanmin(rotarr)):
                vort[h,d]=np.nanmin(rotarr)

    return vort, dist, height


##########
def read_from_radar_file(config, radar_file, is_WSR=False):
    if is_WSR == False:
        radar = pyart.io.read(radar_file)
    elif is_WSR == True:
        if config.radar_controls['Use_downloaded_files']== True:
            radar = pyart.io.read(radar_file)
        else:
            radar = pyart.io.read_nexrad_archive(radar_file)
    return radar
#  function_cache_memory = Memory(plot_config.g_cache_directory,verbose=1)
#  cached_read_from_radar_file = function_cache_memory.cache(read_from_radar_file)
##########
def sweep_index(tilt, is_WSR, radar= None):
    swp_id = None
    if is_WSR == True:
        #Hard code the swp numbers that will be associated with a given tilt angle
        if tilt == 0.5: swp_id=[0 , 1]
        elif tilt == 1.0: swp_id=[2 , 3]
        elif tilt == 1.5: swp_id=[4 , 5]
        else: print('The tilt angle {} is not hard coded yet for WSR'.format(tilt))

    else:
        #  print('radar.get_elevation: '+radar.get_elevation())
        for i in range(radar.nsweeps):
            #  print('Tilt: '+str(tilt))
            #  print('i: '+str(i))

            ## Det the actual tilt angle of a given sweep (returns an array)
            tilt_ang = radar.get_elevation(i)
            #  print('tilt_ang: '+str(tilt_ang))
            ## Check to see if the radarfile matches the elevation tilt we are interested in
            if np.around(tilt_ang[0], decimals=1) == tilt: swp_id = i
    #  print('swp_id: '+str(swp_id))

    return swp_id
##########
def det_nearest_WSR(p_df):
    ''' locate the nearest WSR88D site to the specified insitu instruments
    '''
    #find the locations of all WSR88D sites
    #(outputs dict in format {Site_ID:{lat:...,lon:...,elav:...], ...})
    all_WSR = pyart.io.nexrad_common.NEXRAD_LOCATIONS
    #set up empty dataframe with site Ids as column names
    d_from_all_r = pd.DataFrame(columns = all_WSR.keys())
    #fill in said dataframe with the distance from all 88D sites from each probe measurement
    for key in all_WSR:
        d_from_r = np.square(p_df['lat']-all_WSR[key]['lat']) + np.square(p_df['lon']-all_WSR[key]['lon'])
        d_from_all_r[key] = d_from_r

    #Determine which WS88D site is closest to the probe and add to the original probe dataframe
    p_df['Radar_ID'] = d_from_all_r.idxmin(axis = 1)
    #Determine the unique radar sites to be plotted
    r_ofintrest = p_df.Radar_ID.unique()

    #Set up dict with format{Rsite_ID : [start_oftimerange_for_given_rad : end_oftimerange_for_given_rad]}
    WSR_dict=dict()
    for rad_site in  r_ofintrest:
        trange_r = p_df.loc[p_df.Radar_ID == rad_site, ['datetime']].rename(columns={'datetime': rad_site})
        WSR_dict.update({rad_site:[trange_r.min().get(rad_site), trange_r.max().get(rad_site)]})
    return WSR_dict
##########
def get_WSR_from_AWS(config, day, start, end, radar_id):
    ''' Retrieve the NEXRAD files that fall within a timerange for a specified radar site from the AWS server
    ----------
    INPUTS  radar_id : string,  four letter radar designation
            start & end: datetime,  start & end of the desired timerange
    -------
    RETURN  radar_list : Py-ART Radar Objects
    '''
    # Create this at the point of use # Otherwise it saves everything and eventually crashes
    conn = nexradaws.NexradAwsInterface()
    #Determine the radar scans that fall within the time range for a given radar site
    scans = conn.get_avail_scans_in_range(start, end, radar_id)
    print("There are {} scans available between {} and {}\n".format(len(scans), start, end))

    # Don't download files that you already have...
    path =  config.g_download_directory+ day +'/radar/Nexrad/Nexrad_files/'
    # If you dont have the path already make it and download the files
    if not os.path.exists(path): Path(path).mkdir(parents=True, exist_ok=True)
    # Remove all files ending in _MDM
    scans = list(filter(lambda x: not fnmatch.fnmatch(x.create_filepath(path, False)[-1], '*_MDM') , scans))
    # missing_scans is a list of scans we don't have and need to download
    # create_filepath returns tuple of (directory, directory+filename)
    # [-1] returns the directory+filename
    missing_scans = list(filter(lambda x: not Path(x.create_filepath(path,False)[-1]).exists(), scans))
    # missing files is the list of filenames of files we need to download
    missing_files = list(map(lambda x: x.create_filepath(path,False)[-1], missing_scans))
    print("missing ", len(missing_files), "of ", len(scans), " files")
    print(missing_files)

    # Download the files
    results = conn.download(missing_scans, path, keep_aws_folders=False)
    print(results.success)
    print("{} downloads failed: {}\n".format(results.failed_count,results.failed))
    #print("Results.iter_success : {}\n".format(results.iter_success()))

    # missing_scans_after is a list of scans we don't have (download failed)
    # create_filepath returns tuple of (directory, directory+filename)
    # [-1] returns the directory+filename
    missing_files_after = list(filter(lambda x: not Path(x.create_filepath(path,False)[-1]).exists(), scans))

    if len(missing_files_after) > 0:
        print("ERROR: Some Radar Scans Missing \n ", missing_files_after)
        exit()

    # Return list of files
    radar_files = list(map(lambda x: x.create_filepath(path,False)[-1], scans))
    return radar_files

##########
def smooth_data(radar, data_name):
    """
    Smooth and replace input data using a median filter technique built into scipy
    Parameters:
    ===========
    radar: struct
        pyart radar object
    data_name: string
        name of data field in radar object to be smoothed
    Returns:
    ========
    none
    """
    data           = radar.fields[data_name]['data']
    smooth_data    = scipy.ndimage.filters.median_filter(data, 3)
    smooth_data_ma = np.ma.masked_where(np.ma.getmask(data), smooth_data)
    #  radar.add_field( 'vel_smooth', {'data': smooth_data_ma}, replace_existing = True)
    return smooth_data_ma
def aliasfix(array,delta,nyq):
    '''
   Half-assed attempt to fix de-aliasing errors. Calculates texture
   (difference between point and mean) over a floating window of size delta x delta in the domain.
   If the texture is greater than 1.5*nyquist (arbitrary threshold), then it is a de-aliasing error.
   The if statements fix this according to the nyquist velocity.
    '''

    mean = convolve(array,np.ones((delta,delta))/delta**2.)
    maskpos = np.logical_and(np.abs(array-mean)>nyq*1.5,array>0)
    maskneg = np.logical_and(np.abs(array-mean)>nyq*1.5,array<0)

    array[maskpos]= -2.*nyq+array[maskpos]
    array[maskneg]= 2.*nyq+array[maskneg]

    return array
####
def llsdmain(radar, ref_name, vel_name, swp_id):
    """
    Main processing function for LLSD, applies smoothing and masks before calling llsd compute
    Parameters:
    ===========
    radar: struct
        pyart radar object
    ref_name: string
        name of reflecitivty field
    vel_name: string
        name of doppler velocity field
    Returns:
    ========
    hdr:
        azimuthal shear calculated via the linear least squares derivitives method
    """
    FILLVALUE = -9999
    SCALING = 10

    #define the indices for the required sweep
    sweep_startidx = np.int64(radar.sweep_start_ray_index['data'][swp_id])
    sweep_endidx = np.int64(radar.sweep_end_ray_index['data'][swp_id])

    #  data quality controls on entire volume
    #  smooth_data(radar, vel_nam)

    #extract data
    r= radar.range['data']
    theta = radar.azimuth['data'][sweep_startidx:sweep_endidx+1]
    #  print(theta)
    #  theta = theta*np.pi/180
    #  theta, r = np.meshgrid(theta, r)
    r, theta = np.meshgrid(r, theta)

    refl_full= radar.fields[ref_name]['data']
    vrad_ma  = radar.fields[vel_name]['data'][sweep_startidx:sweep_endidx+1]
    vrad     = np.ma.filled(vrad_ma, fill_value=0)
    mask     = np.ma.getmask(vrad_ma)
    #call llsd compute function
    azi_shear_tilt = lssd_compute(r, theta, vrad, mask, sweep_startidx, sweep_endidx,refl_full,type_grad='az')
    #  azi_shear_tilt = lssd_compute(r, theta, vrad, mask, sweep_startidx, sweep_endidx,refl_full,type_grad='az')
    #scale
    azi_shear_tilt = azi_shear_tilt*SCALING

    #init az_shear for the volume
    azi_shear = np.zeros(refl_full.shape)
    #insert az shear tilt into volume array
    np.set_printoptions(threshold=sys.maxsize)
    azi_shear[sweep_startidx:sweep_endidx+1] = azi_shear_tilt
    #  print(azi_shear)
    #define meta data
    azi_shear_meta = {'data': azi_shear,
                      'long_name': 'LLSD Azimuthal Shear',
                      '_FillValue': FILLVALUE,
                      '_Least_significant_digit': 2,
                      'comment': f'Scaled by x{SCALING}. LLSD azimuthal shear calculation from Miller, M. L., Lakshmanan, V., and Smith, T. M. (2013). An Automated Method for Depicting Mesocyclone Paths and Intensities. Weather and Forecasting, 28(3): 570-585. Effective range of this technique is limited by the window size.',
                      'units': 'Hz'}
    #return shear data
    return azi_shear_meta

def lssd_compute(r_fromradar, theta, vrad, mask, sweep_startidx, sweep_endidx, refl_full, type_grad):
    """
    Compute core for llsd, uses numpy only functions and numba jit.
    Parameters:
    ===========
    r: array
        volume radar range array (2D) (m)
    theta: array
        volume radar azimuth angle array (2D) (radians)
    vrad: array
        volume radial velocity array (2D) (m/s)
    mask: logical array
        volume mask for valid radial velocities (2D)
    sweep_startidx: numba int64 array
        index of starting rays for tilts
    sweep_endidx: numba int64 array
        index of ending rays for tilts

    Returns:
    ========
    azi_shear:
        azimuthal shear calculated via the linear least squaresderivitives method
    """
    #set the constants definining the LLSD grid in the azimuthal and radial directions
    if type_grad == 'az':
        azi_saxis = 500#m              #notes: total azimuthal size = 2*azi_saxis
        rng_saxis = 1  #idx away from i  #notes: total range size = 2*rng_saxis
    elif type_grad == 'div':
        azi_saxis = 500#m              #notes: total azimuthal size = 2*azi_saxis
        rng_saxis = 1  #idx away from i  #notes: total range size = 2*rng_saxis
    elif type_grad =='div_19':
        azi_saxis = 250#m              #notes: total azimuthal size = 2*azi_saxis
        rng_saxis = 3  #idx away from i  #notes: total range size = 2*rng_saxis
    elif type_grad =='az_19':
        azi_saxis = 500#m              #notes: total azimuthal size = 2*azi_saxis
        rng_saxis = 1  #idx away from i  #notes: total range size = 2*rng_saxis

    #get size and init az_shear_tilt
    sz = vrad.shape
    #  print(vrad.shape)
    shear_tilt= np.zeros(sz)

    #  theta = np.ma.filled(theta, fill_value=0)
    #  r_fromradar= np.ma.filled(r_fromradar, fill_value=0)

    #determine the az extent of the [half] kernal along the radial
    az_k_datapnts=[0] #the zero here just shifts the index so its in line with c_rangebin later on (everything is shifted +1)
    for c_rangebin in np.arange(0 + rng_saxis, sz[1] - rng_saxis): #iterate over each range bin
        distfrom_rad_to_cp = r_fromradar[0, c_rangebin]
        center_raytheta= theta[0,c_rangebin]
        max_offset_theta= math.degrees(math.atan(azi_saxis/distfrom_rad_to_cp))

        #number of kernals (for half kernal)
        num_k =0
        for ray in np.arange(0, sz[0]): #iterate over each ray
            ray_theta = theta[ray, 0]
            offset_theta= ray_theta-center_raytheta
            #  print('ray:{} center:{} offset:{} maxoffset:{}'.format(ray_theta,
                                        #  center_raytheta, offset_theta, max_offset_theta))
            if abs(offset_theta) <= abs(max_offset_theta):
                num_k = num_k+1
            else:
               break

        #limit the offset to 100(100 is the legacy)
        if num_k> 100:  num_k= 100

        az_k_datapnts.append(num_k)
    #  print(az_k_datapnts)

    for c_ray in np.arange(0, sz[0]): #iterate over each ray
        for c_rangebin in np.arange(0 + rng_saxis, sz[1] - rng_saxis): #iterate over each range bin
            #skip if centered bin is masked
            if mask[c_ray, c_rangebin]: pass
            else:
                #define the indices for the LLSd grid and deal with wrapping
                #  print('c_ray:{} c_rangebin:{} '.format(c_ray, c_rangebin))
                j_max, j_min = c_ray+abs(az_k_datapnts[c_rangebin]), c_ray-abs(az_k_datapnts[c_rangebin])
                i_max, i_min = c_rangebin+rng_saxis, c_rangebin-rng_saxis
                if j_min< 0: j_min = j_min + sz[0]
                if j_max > sz[0]-1: j_max = j_max - sz[0]
                #  print('j_max:{} j_min:{} '.format(j_max, j_min))
                #  print('i_max:{} i_min:{} '.format(i_max, i_min))

                if j_max< j_min:
                    kernal_thetas=theta
                    kernal_thetas=np.delete(kernal_thetas, np.s_[j_max:j_min], axis=0)

                    kernal_ranges=r_fromradar
                    kernal_ranges= np.delete(kernal_ranges, np.s_[j_max:j_min], axis=0)

                    kernal_vrad=vrad
                    kernal_vrad= np.delete(kernal_vrad, np.s_[j_max:j_min], axis=0)
                else:
                    kernal_thetas=theta[j_min:j_max+1]
                    kernal_ranges=r_fromradar[j_min:j_max+1]
                    kernal_vrad=vrad[j_min:j_max+1]

                #perform calculations according to Miller et al., (2013)
                topsum, botsum = 0, 0
                masked = False
                if type_grad == 'az':
                    dtheta = (kernal_thetas[:,i_min:i_max+1] - theta[c_ray, c_rangebin])
                    topsum = np.sum((kernal_ranges[:, i_min:i_max+1]*dtheta) * kernal_vrad[:, i_min:i_max+1])
                    botsum = np.sum((kernal_ranges[:, i_min:i_max+1]*dtheta)**2)
                    #arclength = 2pi(r)(theta/360)
                    if mask[c_ray, c_rangebin]: masked = True

                elif type_grad == 'div':
                    drange= 15
                    i_range= np.arange(i_min, i_max+1)
                    i =  i_range-c_rangebin
                    topsum = np.sum(i * kernal_vrad[:, i_min:i_max+1])
                    botsum = drange*(np.sum((i)**2))
                elif type_grad in ['div_19', 'az_19']:
                    dr= (kernal_ranges[:, i_min:i_max+1]-r_fromradar[c_ray, c_rangebin])
                    dtheta = (kernal_thetas[:,i_min:i_max+1] - theta[c_ray, c_rangebin])
                    u = kernal_vrad[:, i_min:i_max+1]

                    detA=np.sum(dr**2)*np.sum(dtheta**2)
                    detB=(-1)*(np.sum(dtheta**2)*np.sum(dr)*np.sum(dr))
                    detC=(2)*(np.sum(dr*dtheta)*np.sum(dr)*np.sum(dtheta))
                    detD=(-1)*(np.sum(dr*dtheta)*np.sum(dr*dtheta))
                    det=(detA+detB+detC+detD)**(-1)

                    if type_grad == 'div_19':
                         inner_termA= (np.sum(dr)*np.sum(dtheta))-(np.sum(dr*dtheta))
                         inner_termB= (np.sum(dr**2))-(np.sum(dr)*np.sum(dr))
                         inner_termC= (np.sum(dr*dtheta)*np.sum(dr))-(np.sum(dtheta)*np.sum(dr**2))
                    elif type_grad == 'az_19':
                         inner_termA= (np.sum(dr**2)-(np.sum(dtheta*dtheta)))
                         inner_termB= (np.sum(dr)*np.sum(dtheta))-(np.sum(dr*dtheta))
                         inner_termC= (np.sum(dr*dtheta)*np.sum(dtheta**2))-(np.sum(dr)*np.sum(dtheta**2))
                    termA= (np.sum(dr*u))*inner_termA
                    termB= (np.sum(dtheta*u))*inner_termB
                    termC= (np.sum(u))*inner_termC
                    #  print('sum u:{} sum dt*u:{} sum dr*u:{}'.format(np.sum(u), np.sum(dtheta*u), np.sum(dr*u)))
                    #  print('inA:{} inB:{} inC:{}'.format(inner_termA, inner_termB, inner_termC))
                    #  print('det:{} u:{} dr:{}'.format(det, u, dr))
                    #  print('termA:{} termB:{} termC:{}'.format(termA, termB, termC))
                #exclude regions which contain any masked pixels
                if masked:
                    #  print('masked')
                    pass
                else:
                    #  print('hey')
                    #  print('topsum:{} botsum:{}'.format(topsum, botsum))
                    #  print('top/bot:{}'.format(topsum/botsum))
                    if type_grad in ['az','div']:
                        #exclude areas where there is only one point in each grid
                        if botsum == 0:
                            pass
                        else:
                            shear_tilt[c_ray, c_rangebin] = topsum/botsum
                    elif type_grad in ['az_19', 'div_19']:
                        shear_tilt[c_ray, c_rangebin]= det*(termA+termB+termC)
                #  if c_ray in range(300,400):
                    #  azi_shear_tilt[c_ray, c_rangebin] = 300
        #  print(shear_tilt)
    return shear_tilt


'''
    #begin looping over grid
    for i in np.arange(0, sz[0]): #iterate over each ray
        for j in np.arange(0 + rng_saxis, sz[1] - rng_saxis): #iterate over each range bin
            #skip if j is invalid
            if mask[i, j]:
                continue
            #defining the amount of index offsets for azimuthal direction
            arc_len_idx_offset = np.int64([(azi_saxis//((2*r_fromradar[i, j]*np.pi)/360))])[0] #arc length as a fraction or circ
            #limit the offset to 100
            #  if arc_len_idx_offset > 100:
                #  arc_len_idx_offset = 100
            #define the indices for the LLSd grid and deal with wrapping
            lower_arc_idx = i - arc_len_idx_offset
            upper_arc_idx = i + arc_len_idx_offset
            if lower_arc_idx < 0: lower_arc_idx = lower_arc_idx + sz[0]
            if upper_arc_idx > sz[0]-1: upper_arc_idx = upper_arc_idx - sz[0]
            if upper_arc_idx < lower_arc_idx:
                ii_range = np.concatenate((np.arange(lower_arc_idx, sz[0], 1), np.arange(0, upper_arc_idx+1 ,1)), axis=0)
            else:
                ii_range = np.arange(lower_arc_idx, upper_arc_idx+1)
            #define jj range
            jj_range = np.arange(j-rng_saxis+1, j+rng_saxis)
            #perform calculations according to Miller et al., (2013)
            topsum, botsum = 0, 0
            masked = False
            for ii in ii_range:
                for jj in jj_range:
                    if type_grad == 'az':
                        dtheta = (theta[ii, jj] - theta[i, j])
                        #ensure the angle difference doesnt wrap onto another tilt
                        if (abs(dtheta) > np.pi) and (dtheta > 0):
                            dtheta = ((theta[ii, jj]-2*np.pi) - theta[i, j])
                        elif (abs(dtheta) > np.pi) and (dtheta < 0):
                            dtheta=(theta[ii, jj]) - (theta[i, j]-2*np.pi)
                        topsum = topsum + (r_fromradar[ii, jj]*dtheta) * vrad[ii, jj]
                        botsum = botsum + (r_fromradar[ii, jj]*dtheta)**2
                        if mask[ii, jj]:
                            masked = True
                    elif type_grad == 'div':
                        pass


            if masked:
                #exclude regions which contain any masked pixels
                pass
            elif botsum == 0:
                #exclude areas where there is only one point in each grid
                pass
            else:
                azi_shear_tilt[i, j] = topsum/botsum
            if i in range(300,400):
                azi_shear_tilt[i, j] = 300

    #init az_shear for the volume
    azi_shear = np.zeros(refl_full.shape)
    #insert az shear tilt into volume array
    azi_shear[sweep_startidx:sweep_endidx+1] = azi_shear_tilt

    return azi_shear
'''
##############################################################
###############################
### Command line prompted Defns
###############################
def info_radar_files(level, files):
    for radar_file in files[:]:
        print("\nInfo ... " + " Input File: " + str(radar_file))
        radar = pyart.io.read(radar_file)
        radar.info(level = level)
        print(radar.info(level='compact'))
        #  print('Nyquist Vel: {}'.format(radar.get_nyquist_vel(swp_id)))
        #  print('Vel limits: {}'.format(pyart.config.get_field_limits('velocity', container=radar)))
        #  print('Refl limits: {}'.format(pyart.config.get_field_limits('reflectivity', container=radar)))
        #  print('range: {}'.format(radar.range))
        #  print('Gate Area: {}'.format(radar.get_gate_area(swp_id)))
        #  print('XYZ: {}'.format(radar.get_gate_x_y_z(swp_id)))


#*****
def dealias_radar(radar, vel_field, texture_threshold, is_WSR=False):
    #  texture_feild= pyart.retrieve.calculate_velocity_texture(radar, vel_field=vel_field, wind_size=3)
    #  radar.add_field('vel_texture', texture_feild, replace_existing=True)

    if  texture_threshold != False:
        texture_feild= pyart.retrieve.calculate_velocity_texture(radar, vel_field=vel_field, wind_size=3)
        radar.add_field('vel_texture', texture_feild, replace_existing=True)

        gatefilter = pyart.filters.GateFilter(radar)
        #  ncp_name, ncp_thresh = 'normalized_coherent_power', .4
        #  gatefilter.exclude_below(ncp_name, ncp_thresh)
        gatefilter.exclude_above('vel_texture', texture_threshold)
        #  gatefilter.exclude_equal('vel_texture', 0)
        if is_WSR== False:
            dealias_data = pyart.correct.dealias_region_based(radar, vel_field=vel_field, gatefilter=gatefilter)
    else:
        if is_WSR== False:
            dealias_data = pyart.correct.dealias_region_based(radar, vel_field=vel_field)
        elif is_WSR == True:
            dealias_data = pyart.correct.dealias_unwrap_phase(radar, vel_field=vel_field)#, unwrap_unit='volume')
    return dealias_data, radar

def dealias_radar_file(prefix, radar_file, vel_field, texture_threshold, is_WSR=False):
    out_filename = os.path.dirname(radar_file) + '/' + prefix + os.path.basename(radar_file)
    print("\nDealiasing... " + " Input File: " + str(radar_file) + " Output File: " + out_filename)

    if radar_file == '.':
        return;

    radar = pyart.io.read(radar_file)

    if radar.scan_type == 'rhi':
        print('Skipping, we are not dealiasing RHI files at this time \n')
        return;

    dealias_data, radar = dealias_radar( radar, vel_field, texture_threshold )
    radar.add_field('corrected_velocity', dealias_data, replace_existing=True)

    pyart.io.write_cfradial(out_filename, radar)


def dealias_radar_files(prefix, files, Rtype):
    if Rtype == 'KA':
        vel_field = 'velocity'
        texture_thresh= False
        is_WSR=False
    elif Rtype == 'NOXP':
        vel_field = 'VEL'
        texture_thresh= 3
        is_WSR=False
    elif Rtype == 'WSR':
        vel_field = 'velocity'
        texture_thresh= False
        is_WSR=True

    #for radar_file in files[:]:
    Parallel(n_jobs=-1, verbose=10) (delayed(dealias_radar_file)(prefix, radar_file, vel_field, texture_thresh, is_WSR) for radar_file in files[:])

def fix_NOXP_order(prefix, files):
    ''' pyart requires rays in sweep to be ordered by angle
        Some radars sweeps cross the 360 degree threshold (359->0->1)
        Some radars sweeps contain segments that are out of order
        We need to fix both problems for pyart so all rays in sweep are ordered by ascending azimuth angle
        *** both azimuth angle array and field arrays are reordered in rfile
        ---
        input:  azimuth_angles = [..., 358, 359,   0,   1,   2, ...,  39,  40, 200, 201, ...]
        output: azimuth_angles = [200, 201, ..., 358, 359, 360, 361, 362, ..., 399, 400, ...]
    '''
    def remove_zero_crossing(azimuth_angles):
        ''' Use angles >= 360 when sweep passes from 359.x to 0.x degrees: Dont use offset on other segments
            ---
            input:  azimuth_angles = [..., 358, 359,   0,   1,   2, ...,  39,  40, 200, 201, ...]
            output: azimuth_angles = [..., 358, 359, 360, 361, 362, ..., 399, 400, 200, 201, ...]
            Note the 359->0 transition and the 40->200 transition
        '''
        previous_angle, zero_crossed = 0, False
        for index, angle in enumerate(azimuth_angles):
            if previous_angle > 350 and angle < 1: zero_crossed = True
            discontinuity = (angle - previous_angle) > 10

            if zero_crossed and discontinuity: zero_crossed = False
            previous_angle = angle

            if zero_crossed: a = angle + 360.0
            else: a = angle
            azimuth_angles[index] = a
        return azimuth_angles

    def sort_on_ascending_angles(azimuth_angles):
        ''' Generate array of indexes into original array that can be used to generate new sorted array
            ---
            input: azimuth_angles = [..., 358, 359, 360, 361, 362, ..., 399, 400, 200, 201, ...]
            output: sorte_index_list = [ 41, 42, ..., 0, 1, 2, ...]
        '''
         #Output:{0: 0, 1: .5, 2: 1, ... n: 359.5}
        index_angle_dict = dict(enumerate(azimuth_angles))
        # Create list of tuples (index, angle) sorted by angle
        sorted_index_angle_tuples = sorted(index_angle_dict.items(), key=lambda kv: kv[1])
        # source index list - Used to reconstruct a list sorted by angles
        sorted_index_list = [a_tuple[0] for a_tuple in sorted_index_angle_tuples]
        # sorted list of all ray angles
        #sorted_angle_list = [a_tuple[1] for a_tuple in sorted_index_angle_tuples]
        return sorted_index_list

    def reorder_array(array, sorted_index_list):
        ''' Reorder array as directed by sorted_index_list
            ---
            input:  array = [4, 5, 6, 1, 2, 3]
                    sorted_index_list = [3, 4, 5, 0, 1, 2]
            output: array = [1, 2, 3, 4, 5, 6]
        '''
        array_copy = copy.deepcopy(array)
        for target_index, source_index in enumerate(sorted_index_list):
          array[target_index] = array_copy[source_index]
        return array

    ### read radar file
    # # # # # # # # # #
    for radar_file in files[:]:
        print(radar_file)
        rfile= pyart.io.read(radar_file)

        ### order_rays_by_angle
        # # # #  # # # #  # # #
        for sweep_id in range(rfile.nsweeps):
            azimuth_angles = rfile.get_azimuth(sweep_id)
            azimuth_angles = remove_zero_crossing(azimuth_angles)
            sorted_index_list = sort_on_ascending_angles(azimuth_angles)

            # reorder azimuth angles so rays are in order of assending angles
            azimuth_angles = reorder_array(azimuth_angles, sorted_index_list)
            for field_name in list(rfile.fields):
                # Ray data array were overwriting
                field_sweep_data_reference = rfile.get_field(sweep_id, field_name, copy=False)
                # reorder field ray data so rays are in order of assending angles
                field_sweep_data_reference = reorder_array(field_sweep_data_reference, sorted_index_list)


        ### write out the corrected radar file
        # # # # # # # # # # # # # # # # # # # #
        out_filename = os.path.dirname(radar_file) + '/' + prefix + os.path.basename(radar_file)
        print("\nSorting... " + " Input File: " + str(radar_file) + " Output File: " + out_filename)
        pyart.io.write_cfradial(out_filename, rfile)

############################
class Radar_File_Functions(object):
    '''inputs from command line for the various code functionalities; thanks to the fire command below
    (aka this code does a number of different things depending what inputs you give it)'''
    def dealias_files(self, *input_files, Rtype, output_file_prefix: str = 'dealiased_'):
        ''' Dealias radar file(s).
            ---
            input_files : array, Radar files to process.
            Rtype: str (KA or NOXP), Type of radar to process (to determine velocity field name)
            output_file_prefix : string, Path prefix to start all output files with.
        '''
        dealias_radar_files(output_file_prefix, input_files, Rtype)

    def dump_file_info(self, *input_files, pyart_info_level: str = 'standard'):
        ''' Describe radar file(s).
            ---
            input_files : array, Radar files to process.
            pyart_info_level : string, Pyart specification for info level.
        '''
        info_radar_files(pyart_info_level, input_files)
    def fix_NOXP_files(self, *input_files, output_file_prefix: str = 'ordered_'):
        ''' Fix the ordering of the rays in the NOXP radar file(s).
            ---
            input_files : array, Radar files to process.
            pyart_info_level : string, Pyart specification for info level.
        '''
        fix_NOXP_order(output_file_prefix, input_files)

############################
def main():
  fire.Fire(Radar_File_Functions)
if __name__ == '__main__':
    main()

