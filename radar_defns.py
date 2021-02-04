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
