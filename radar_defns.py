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
from os.path import expanduser
from pathlib import Path
from joblib import Memory, Parallel, delayed
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
pp = pprint.PrettyPrinter(indent=4)


def radar_fields_prep(config, rfile, radar_type):
    if radar_type == 'KA':
        vel_name, refl_name = 'corrected_velocity', 'reflectivity'
        ncp_name, ncp_thresh = 'normalized_coherent_power', .4
        rhv_name= 'None'

    if radar_type == 'NOXP':
        vel_name , refl_name = 'VEL', 'DBZ'
        ncp_name, ncp_thresh = 'SQI', .5 
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
        if m == 'sim_vel':
            sd_test = singledop.SingleDoppler2D(L=30.0, radar=rfile, range_limits=[0, 20],
                                    sweep_number=0, name_vr='velocity', thin_factor=[4, 12])
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
            vel_field = rfile.fields['vel_fix']['data']
            print(type(vel_field))
            print(np.shape(vel_field))
            vel_grad= np.gradient(vel_field, axis =1)
            print(type(vel_grad))
            rfile.add_field('vel_gradient', {'data': vel_grad}, replace_existing=True) 
            #  rfile=mask(rfile, gatefilter, vel_grad, 'vel_gradient')

    return rfile


##########
def order_rays_by_angle(sweep_id, rfile):
    ''' pyart requires rays in sweep to be ordered by angle
        Some radars sweeps cross the 360 degree threshold (359->0->1)
        Some radars sweeps contain segments that are out of order
        We need to fix both problems for pyart so all rays in sweep are ordered by ascending azimuth angle
        ---
        input:  azimuth_angles = [..., 358, 359,   0,   1,   2, ...,  39,  40, 200, 201, ...]
        output: azimuth_angles = [200, 201, ..., 358, 359, 360, 361, 362, ..., 399, 400, ...]

        both azimuth angle array and field arrays are reordered in rfile
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

def texture(radar, var):
    """ Determine a texture field using an 11pt stdev
    texarray=texture(pyradarobj, field). """
    fld = radar.fields[var]['data']
    print(fld.shape)
    tex = np.ma.zeros(fld.shape)
    for timestep in range(tex.shape[0]):
        ray = np.ma.std(rolling_window(fld[timestep, :], 11), 1)
        tex[timestep, 5:-5] = ray
        tex[timestep, 0:4] = np.ones(4) * ray[0]
        tex[timestep, -5:] = np.ones(5) * ray[-1]
    return tex

def texture_along_ray(radar, var, wind_size=7):
    """
    Compute field texture along ray using a user specified
    window size.

    Parameters
    ----------
    radar : radar object
        The radar object where the field is.
    var : str
        Name of the field which texture has to be computed.
    wind_size : int, optional
        Optional. Size of the rolling window used.

    Returns
    -------
    tex : radar field
        The texture of the specified field.

    """
    half_wind = int((wind_size-1)/2)
    fld = radar.fields[var]['data']
    tex = np.ma.zeros(fld.shape)
    for timestep in range(tex.shape[0]):
        ray = np.ma.std(rolling_window(fld[timestep, :], wind_size), 1)
        tex[timestep, half_wind:-half_wind] = ray
        tex[timestep, 0:half_wind] = np.ones(half_wind) * ray[0]
        tex[timestep, -half_wind:] = np.ones(half_wind) * ray[-1]
    return tex
def angular_texture_2d(image, N, interval):
    """
    Compute the angular texture of an image. Uses convolutions
    in order to speed up texture calculation by a factor of ~50
    compared to using ndimage.generic_filter.

    Parameters
    ----------
    image : 2D array of floats
        The array containing the velocities in which to calculate
        texture from.
    N : int
        This is the window size for calculating texture. The texture will be
        calculated from an N by N window centered around the gate.
    interval : float
        The absolute value of the maximum velocity. In conversion to
        radial coordinates, pi will be defined to be interval
        and -pi will be -interval. It is recommended that interval be
        set to the Nyquist velocity.

    Returns
    -------
    std_dev : float array
        Texture of the radial velocity field.

    """
    # transform distribution from original interval to [-pi, pi]
    interval_max = interval
    interval_min = -interval
    half_width = (interval_max - interval_min) / 2.
    center = interval_min + half_width

    # Calculate parameters needed for angular std. dev
    im = (np.asarray(image) - center) / (half_width) * np.pi
    x = np.cos(im)
    y = np.sin(im)

    # Calculate convolution
    kernel = np.ones((N, N))
    xs = signal.convolve2d(x, kernel, mode="same", boundary="symm")
    ys = signal.convolve2d(y, kernel, mode="same", boundary="symm")
    ns = N**2

    # Calculate norm over specified window
    xmean = xs/ns
    ymean = ys/ns
    norm = np.sqrt(xmean**2 + ymean**2)
    std_dev = np.sqrt(-2 * np.log(norm)) * (half_width) / np.pi
    return std_dev
def calculate_velocity_texture(radar, vel_field=None, wind_size=4, nyq=None,
                               check_nyq_uniform=True):
    """
    Derive the texture of the velocity field.

    Parameters
    ----------
    radar: Radar
        Radar object from which velocity texture field will be made.
    vel_field : str, optional
        Name of the velocity field. A value of None will force Py-ART to
        automatically determine the name of the velocity field.
    wind_size : int, optional
        The size of the window to calculate texture from. The window is
        defined to be a square of size wind_size by wind_size.
    nyq : float, optional
        The nyquist velocity of the radar. A value of None will force Py-ART
        to try and determine this automatically.
    check_nyquist_uniform : bool, optional
        True to check if the Nyquist velocities are uniform for all rays
        within a sweep, False will skip this check. This parameter is ignored
        when the nyq parameter is not None.

    Returns
    -------
    vel_dict: dict
        A dictionary containing the field entries for the radial velocity
        texture.

    """
    # Parse names of velocity field
    if vel_field is None:
        vel_field = get_field_name('velocity')

    # Allocate memory for texture field
    vel_texture = np.zeros(radar.fields[vel_field]['data'].shape)

    # If an array of nyquist velocities is derived, use different
    # nyquist velocites for each sweep in texture calculation according to
    # the nyquist velocity in each sweep.

    if nyq is None:
        # Find nyquist velocity if not specified
        nyq = [radar.get_nyquist_vel(i, check_nyq_uniform) for i in
               range(radar.nsweeps)]
        for i in range(0, radar.nsweeps):
            start_ray, end_ray = radar.get_start_end(i)
            inds = range(start_ray, end_ray)
            vel_sweep = radar.fields[vel_field]['data'][inds]
            vel_texture[inds] = angular_texture_2d(
                vel_sweep, wind_size, nyq[i])
    else:
        vel_texture = angular_texture_2d(
            radar.fields[vel_field]['data'], wind_size, nyq)
    vel_texture_field = get_metadata('velocity')
    vel_texture_field['long_name'] = 'Doppler velocity texture'
    vel_texture_field['standard_name'] = (
        'texture_of_radial_velocity' + '_of_scatters_away_from_instrument')
    vel_texture_field['data'] = ndimage.filters.median_filter(vel_texture,
                                                              size=(wind_size,
                                                                    wind_size))
    return vel_texture_field

