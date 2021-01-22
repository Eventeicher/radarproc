#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.patheffects as PE
from matplotlib.collections import PatchCollection
from matplotlib import ticker
from matplotlib.ticker import (LinearLocator, FixedLocator, MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import ListedColormap, Normalize
from cycler import cycler
from mpl_toolkits.axes_grid1.axes_divider import make_axes_area_auto_adjustable
import matplotlib.cm as cmx
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import LineCollection
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
import pandas as pd
from pandas.plotting import register_matplotlib_converters
import numpy as np
import xarray as xr
from netCDF4 import num2date
import argparse, cProfile, logging, time, os, os.path
from os.path import expanduser
from pathlib import Path
from joblib import Memory, Parallel, delayed
from scipy import ndimage, interpolate
from operator import attrgetter
from collections import namedtuple
import pyart, nexradaws, sys, traceback, shutil, glob, gc, cmocean
import fnmatch
import gc
from pympler.asizeof import asizeof
#this is the file with the plotting controls to access any of the vars in that file use config.var
import config as plot_config
import copy

from read_pforms import read_TInsitu, pform_names
def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))

pd.set_option('display.max_columns', None)
df_holder, p_holder =[], [] 
for p in pform_names('TInsitu'):
    print(p)
    print('............')
    d_avail=read_TInsitu(plot_config, p, True)
    if d_avail == True:
        df, ptype =read_TInsitu(plot_config, p)
        df['pname'] = p
        df = df.drop(columns=['U','V', 'dir', 'spd' ])
        if ptype =='NSSL':
            df = df.drop(columns=['id', 'qc1', 'qc2', 'qc3', 'qc4', 'all_qc_flags'])
        df_holder.append(df)
        p_holder.append(p)
full_df= pd.concat(df_holder)
new_cols=['dist_Prb1', 'dist_Prb2', 'dist_FFld', 'dist_Wins', 'dist_LiDR', 'dist_CoMeT1', 'dist_CoMeT2']
for c in new_cols:
    full_df.insert(2, c, np.nan)
full_df.set_index(['pname', 'datetime'], inplace=True)
#  print(full_df.loc['FFld','datetime'])

test_dict= {}
for pform in p_holder:
    p_df = full_df.loc[pform].index
    for p_time in p_df:
        #actual timeboundary
        time_offset_outer= timedelta(minutes=5)
        #  start_t_outer, end_t_outer = (p_time-time_offset_outer).to_pydatetime(), (p_time + time_offset_outer).to_pydatetime()
        start_t_outer, end_t_outer = (p_time-time_offset_outer), (p_time + time_offset_outer)
        #most data points will fall within the time bound this tighter bound just decreases the #of calcs needed
        time_offset_inner= timedelta(minutes=1)
        #  start_t_inner, end_t_inner = (p_time-time_offset_inner).to_pydatetime(), (p_time + time_offset_inner).to_pydatetime()
        start_t_inner, end_t_inner = (p_time-time_offset_inner), (p_time + time_offset_inner)

        for comp_pform in p_holder:
            #subset data to a smaller time range 
            comp_df = full_df.loc[comp_pform].loc[start_t_inner:end_t_inner]
            if len(comp_df) == 0: 
                #  test if there is data within the actual time window we want (might just need to include more calcs)
                comp_df = full_df.loc[comp_pform].loc[start_t_outer:end_t_outer]
            
            #if there still is not data in the range pass otherwise add to the dataframe
            if len(comp_df) != 0: 
                comp_p_time = nearest(comp_df.index, p_time)
                if comp_pform == 'FFld': comp_col= 'diff_FFld'
                elif comp_pform == 'WinS': comp_col= 'diff_Wins'
                elif comp_pform == 'LIDR': comp_col= 'diff_LiDR'
                elif comp_pform == 'Prb1': comp_col= 'diff_Prb1'
                elif comp_pform == 'Prb2': comp_col= 'diff_Prb2'
                elif comp_pform == 'CoMeT1': comp_col= 'diff_CoMeT1'
                elif comp_pform == 'CoMeT2': comp_col= 'diff_CoMeT2'
                elif comp_pform == 'CoMeT3': comp_col= 'diff_CoMeT3'

                dist = np.square(comp_df['lat']-full_df.loc[(pform, p_time),'lat'])+ np.square(comp_df['lon']-full_df.loc[(pform, p_time),'lon'])
                full_df[(pform,p_time), comp_col] = dist
                #  test_dict.update(comp_pform=comp_p_time)
            #  else:
                #  test_dict.update(comp_pform= np.nan)
    #  print(test_dict)
    print('hey')
    print(full_df)
        #  for comp_pform in p_holder:
            #  print(p_df)
        #  comp_p_time = nearest(comps.loc[comp_pform].index, p_time)
    #  for p_time in p_df:
        #  for comp_pform in p_holder:
            #  print(p_time)
            #  print(comps)
            #  comp_p_time = nearest(comps.loc[comp_pform].index, p_time)
            #  comp_p_time = nearest(full_df.loc[comp_pform].index, p_time)
            #  300 sec = 5min
            #  if abs(comp_p_time - p_time).seconds <= 300:
                #  test_dict.update({comp_pform: comp_p_time})
        #  print(test_dict)
                #  print("p: {}, c_p: {}, p_time: {}, comp_p_time: {}".format(pform, comp_pform, p_time, comp_p_time))
#  print(test.iloc[3])
#  print(near)
