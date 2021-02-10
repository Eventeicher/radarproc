#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import needed modules
######################
import matplotlib
matplotlib.use('agg')
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
import datetime as dt
from datetime import datetime, date, timedelta
import cartopy
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
import argparse, cProfile, logging, time, os, os.path
import os
from os.path import expanduser
import pathlib
from pathlib import Path
from joblib import Memory, Parallel, delayed
import scipy
from scipy import ndimage, interpolate
from operator import attrgetter
from collections import namedtuple
import pyart, nexradaws, sys, traceback, shutil, glob, gc, cmocean
import cftime # for num2pydate
import pprint
import math
import copy
import pytda
import singledop
import pickle
import fire

pp = pprint.PrettyPrinter(indent=4)


def radar_fields_prep(config, rfile, radar_type, sweep_id):
    if radar_type == 'KA':
        vel_name, refl_name = 'corrected_velocity', 'reflectivity'
        ncp_name, ncp_thresh = 'normalized_coherent_power', .4
        rhv_name= 'None'

    if radar_type == 'NOXP':
        vel_name , refl_name = 'corrected_velocity', 'DBZ'
        #  vel_name , refl_name = 'VEL', 'DBZ'
        ncp_name, ncp_thresh = 'SQI', .6 
        rhv_name = 'RHOHV'

    #creating the mask for attenuation
        #  range_mask = np.zeros(np.shape(reflectivity))
        #  for i in range(0, len(range_mask[:,0])): range_mask[i,:] = rfile.range['data'] > (rfile.range['data'][-1]-1000.)
        #  range_mask = range_mask.astype(bool)
        #  total_mask = [any(t) for t in zip(range_mask.flatten(), normal_mask.flatten())]
        #  gatefilter = pyart.filters.GateFilter(rfile)
        #  gatefilter.exclude_below('DBZ', 10)
    gatefilter= pyart.filters.moment_based_gate_filter(rfile, ncp_field= ncp_name, rhv_field=rhv_name, 
                                                      refl_field=refl_name, min_ncp= ncp_thresh)
    # apply the mask to the necessary feilds
    for m in config.r_mom:
        print(m)
        def textures(rfile, call, orig_field='vel_fix'):
            if call== 'calc_vel':
                texture_feild_almost= pyart.retrieve.calculate_velocity_texture(rfile, vel_field=orig_field, wind_size=3)
                texture_feild = texture_feild_almost['data']
            if call== 'text':
                texture_feild = pyart.util.texture(rfile, 'vel_fix')
            if call== 'aray_text':
                texture_feild = pyart.util.texture_along_ray(rfile, 'vel_fix')
            return texture_feild

        def mask(rfile, gatefilter, moment, new_mom_name):
            mom_masked= np.ma.masked_where(gatefilter.gate_included == False, moment)
            rfile.add_field(new_mom_name, {'data': mom_masked}, replace_existing=True) 
            return rfile

        ###
        if m in ['vel', 'refl']:
            if m == 'vel': orig_field_name, new_feild_name = vel_name, 'vel_fix'
            if m == 'refl': orig_field_name, new_feild_name = refl_name, 'refl_fix'
            orig_field = rfile.fields[orig_field_name]['data']
            rfile=mask(rfile, gatefilter, orig_field, new_feild_name)
        
        if m in ['vel_texture', 'vel_text_texture', 'vel_aray_texture']: 
            if m == 'vel_texture':
                call = 'calc_vel'
            if m == 'vel_text_texture':
                call = 'text'
            if m == 'vel_aray_texture':
                call = 'aray_text'
            texture_feild=textures(rfile, call)
            rfile=mask(rfile, gatefilter, texture_feild, m)
        #  if m == 'sim_vel':
            #  sd_test = singledop.SingleDoppler2D(L=30.0, radar=rfile, range_limits=[0, 20],
                                    #  sweep_number=0, name_vr='velocity', thin_factor=[4, 12])
#
            #  sim_V= pyart.util.simulated_vel_from_profile(rfile, rfile.fields['vel_fix']['data'])
            #  rfile=mask(rfile, gatefilter, sd_test, m)
        if m == 'diff_dbz':
            differnce_feild_name= 'difference'
            tot_field = rfile.fields['DBZ_TOT']['data']
            dbz_field = rfile.fields['DBZ']['data']
            diff = tot_field - dbz_field 
            rfile=mask(rfile, gatefilter, diff, 'difference')
        if m == 'vel_grad':
            #  vel_field = rfile.get_field(sweep_id, 'vel_fix', copy=False)
            vel_field = rfile.fields['vel_fix']['data']
            print(np.shape(rfile.get_field(sweep_id, 'corrected_velocity')))
            #  vel_smoothed=aliasfix(rfile.get_field(sweep_id,'velocity'),13,rfile.get_nyquist_vel(0))
            vel_smoothed=aliasfix(rfile.get_field(sweep_id,'corrected_velocity'),13,rfile.get_nyquist_vel(0))
            #  vel_smoothed = scipy.signal.savgol_filter(vel_field, 15, 1, axis=1)
            vel_grad= np.gradient(vel_smoothed, 15, axis= 1)*100
            #  rfile.add_field('vel_gradient', {'data': vel_grad}, replace_existing=True)
            rfile=mask(rfile, gatefilter, vel_grad, 'vel_gradient')
        if m == 'vel_grad0_0':
            vel_field = rfile.fields['vel_fix']['data']
            vel_grad = scipy.signal.savgol_filter(vel_field, window_length=19, polyorder=2, deriv=0, axis=0)
            rfile=mask(rfile, gatefilter, vel_grad, 'vel_gradient0_0')

        if m == 'vel_grad1_0':
            vel_field = rfile.fields['vel_fix']['data']
            vel_grad = scipy.signal.savgol_filter(vel_field, window_length=19, polyorder=2, deriv=1, axis=0)
            rfile=mask(rfile, gatefilter, vel_grad, 'vel_gradient1_0')

        if m == 'vel_grad2_0':
            vel_field = rfile.fields['vel_fix']['data']
            vel_grad = scipy.signal.savgol_filter(vel_field, window_length=19, polyorder=2, deriv=2, axis=0)
            rfile=mask(rfile, gatefiter, vel_grad, 'vel_gradient2_0')


        if m == 'vel_grad0_1':
            vel_field = rfile.fields['vel_fix']['data']
            vel_grad = scipy.signal.savgol_filter(vel_field, window_length=19, polyorder=2, deriv=0, axis=1)

            rfile=mask(rfile, gatefilter, vel_grad, 'vel_gradient0_1')

        if m == 'vel_grad1_1':
            vel_field = rfile.fields['vel_fix']['data']
            vel_grad = scipy.signal.savgol_filter(vel_field, window_length=19, polyorder=2, deriv=1, axis=1)
            rfile=mask(rfile, gatefilter, vel_grad, 'vel_gradient1_1')

        if m == 'vel_grad2_1':
            vel_field = rfile.fields['vel_fix']['data']
            vel_grad = scipy.signal.savgol_filter(vel_field, window_length=19, polyorder=2, deriv=2, axis=1)
            rfile=mask(rfile, gatefilter, vel_grad, 'vel_gradient2_1')
    return rfile

##########
def read_from_radar_file(radar_file, is_WSR=False):
    if is_WSR == False: 
        radar = pyart.io.read(radar_file)
    elif is_WSR == True: 
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
        for i in range(radar.nsweeps):
            ## Det the actual tilt angle of a given sweep (returns an array)
            tilt_ang = radar.get_elevation(i)
            ## Check to see if the radarfile matches the elevation tilt we are interested in
            if np.around(tilt_ang[0], decimals=1) == tilt: swp_id = i
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
def aliasfix(array,delta,nyq):
    '''
   Half-assed attempt to fix de-aliasing errors. Calculates texture
   (difference between point and mean) over a floating window of size delta x delta in the domain.
   If the texture is greater than 1.5*nyquist (arbitrary threshold), then it is a de-aliasing error.
   The if statements fix this according to the nyquist velocity.
    '''
 
    mean = scipy.convolve(array,np.ones((delta,delta))/delta**2.)
    maskpos = np.logical_and(np.abs(array-mean)>nyq*1.5,array>0)
    maskneg = np.logical_and(np.abs(array-mean)>nyq*1.5,array<0)
 
    array[maskpos]= -2.*nyq+array[maskpos]
    array[maskneg]= 2.*nyq+array[maskneg]
 
    return array
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


def dealias_radar_files(prefix, files, Rtype):
    dealiased_files, duplicated_files, skipped_files = [], [], []
    ###Set Up Temporary Directory
    cachedir='./cachedir'
    mem= Memory(cachedir,verbose=1)
    def dealias_radar(radar_pickled):
        # Regenerate the calculated values in the object
        radar = pickle.loads(radar_pickled)
        dealias_data = pyart.correct.region_dealias.dealias_region_based(radar, vel_field=vel_field)
        return dealias_data

    ###joblib dealias stuff
    #Don't use cached results
    #  persistant_dealias_radar = dealias_radar
    #Use caching
    persistant_dealias_radar = mem.cache( dealias_radar )

    if Rtype == 'KA':
        vel_field = 'velocity'
    elif Rtype == 'NOXP':
        vel_field = 'VEL'

    for radar_file in files[:]:
        #  print("a: {},\n b:{},\n c:{},\n d:{},\n e:{}".format(radar_file, str(radar_file),
                                #  os.path.dirname(radar_file), prefix, os.path.basename(radar_file)))
        out_filename = os.path.dirname(radar_file) + '/' + prefix + os.path.basename(radar_file)
        print("\nDealiasing... " + " Input File: " + str(radar_file) + " Output File: " + out_filename)

        if radar_file == '.':
            skipped_files.append(radar_file)
            continue;

        radar = pyart.io.read(radar_file)

        if radar.scan_type == 'rhi':
            print('Skipping, we are not dealiasing RHI files at this time \n')
            skipped_files.append(radar_file)
            continue;

        #  if 'corrected_velocity' in radar.fields:
            #  print('Input file already contains dealias data in field corrected_velocity......'
                      #  ' so just creating a copy with the specified prefix \n')
            #  duplicated_files.append(out_filename)
        #  else:
        dealias_data = persistant_dealias_radar( pickle.dumps(radar) )
        radar.add_field('corrected_velocity', dealias_data, replace_existing=True)
        dealiased_files.append(out_filename)

        pyart.io.write_cfradial(out_filename, radar)

    print("\n dealiased_files:")
    pprint(dealiased_files)
    print("\n skipped_files:")
    pprint(skipped_files)
    print("\n duplicated_files:")
    pprint(duplicated_files)


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

