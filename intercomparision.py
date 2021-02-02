#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib
import tabulate
import mpu
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
pd.set_option('display.max_columns', None)

## Whether to calc the distance between pforms 
make_dist_csv= False
 #  True
apply_dist_threshold = False
apply_storm_threshold = True#False
plotting = True
comp_plots, prox_plots = True, True
#True#False
#  False
# # # # # # #

## Thresholds
# # # # # # #
#frequency at which to bin (average) the obs
bined_data_freq = '10S'
#threshold between two obs for valid comparisions
valid_time_diff= 5 #Min
valid_dist_diff= 10 #'10km'
#directory path 
#  day_list= ['20190517', '20190518','20190520', '20190523', '20190524','20190525','20190526', '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
day_list= ['20190517', '20190518', '20190523', '20190524','20190525','20190526', '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
day_cuttoff= {'20190517':'23:00', '20190518':'22:30','20190520':'17:00', '20190523':'22:00', '20190524':'19:00','20190525':'19:00','20190526':'18:00',
              '20190527':'20:00', '20190528':'22:15', '20190608':'20:30', '20190611':'22:50', '20190613':'00:00', '20190615':'22:00'}
colors ={'Prb1': 'red', 'Prb2':'green', 'WinS':'blue', 'LIDR': 'black', 'FFld':'gray', 'CoMeT1':'purple', 'CoMeT2':'pink', 'CoMeT3':'orange'}
#  day_list=['20190523']
for day in day_list:
    outdir = plot_config.g_TORUS_directory+day+'/data/mesonets/'

    ## Definitions 
    # # # # # # #
    def nearest(items, pivot):
        return min(items, key=lambda x: abs(x - pivot))

    ## Make a csv with the dist pforms are from each other
    # # # # #  # # # # # # # # # # # # # # # # # # # # # #
    if make_dist_csv == True:
        df_holder, p_holder =[], [] 
        #  pform_list=["CoMeT1", "Prb2" ]
        for p in pform_names('TInsitu'):
            d_avail=read_TInsitu(plot_config, day, p, True)
            if d_avail == True:
                #read in data
                df = pd.DataFrame.empty
                df, ptype =read_TInsitu(plot_config, day, p)
                
                #remove unneeded data columns
                df = df.drop(columns=['U','V', 'dir', 'spd' ])
                if ptype =='NSSL':
                    #  mask_allqc_df = np.ma.masked_where(df['all_qc_flags'].values > 0, df[var].values) #add masked dataset to the object
                    #  mask_allqc_df = np.ma.masked_where(df['qc1'].values > 0, df[var].values) #add masked dataset to the object
                    #  mask_allqc_df = np.ma.masked_where(df['qc4'].values > 0, df[var].values) #add masked dataset to the object
                    #  df = df[df['all_qc_flags']==0]
                    df = df[df['qc1']==0]
                    df = df[df['qc4']==0]
                    df = df.drop(columns=['id', 'qc1', 'qc2', 'qc3', 'qc4', 'all_qc_flags'])
                
                #group the obs together in x second intervals 
                df=df.groupby(pd.Grouper(key="datetime", freq = bined_data_freq)).mean()
                df.reset_index(inplace=True)
                #drop rows with no data 
                df.dropna(subset=['lat', 'lon'], inplace=True)
                #add column with the instrument name
                df.loc[:,'pname'] = p

                #add to the list to be concationated
                df_holder.append(df)
                p_holder.append(p)
        full_df= pd.concat(df_holder)
        full_df.set_index(['pname', 'datetime'], inplace=True)

        # * * * * * *
        ## calc distance between platforms for each data entry
        for pform in p_holder:
            p_df = full_df.loc[pform].index
            for p_time in p_df:
                #actual timeboundary
                time_offset_outer= timedelta(minutes = valid_time_diff)
                #  start_t_outer, end_t_outer = (p_time-time_offset_outer).to_pydatetime(), (p_time + time_offset_outer).to_pydatetime()
                start_t_outer, end_t_outer = (p_time-time_offset_outer), (p_time + time_offset_outer)
                #most data points will fall within the time bound this tighter bound just decreases the #of calcs needed
                time_offset_inner= timedelta(minutes=1)
                #  star t_t_inner, end_t_inner = (p_time-time_offset_inner).to_pydatetime(), (p_time + time_offset_inner).to_pydatetime()
                start_t_inner, end_t_inner = (p_time-time_offset_inner), (p_time + time_offset_inner)

                for comp_pform in p_holder:
                    ##subset data to a smaller time range 
                    comp_df = full_df.loc[comp_pform].loc[start_t_inner:end_t_inner]
                    if len(comp_df) == 0: 
                        #  test if there is data within the actual time window we want (might just need to include more calcs)
                        comp_df = full_df.loc[comp_pform].loc[start_t_outer:end_t_outer]
                    # ***

                    ##if there still is not data in the range pass otherwise add to the dataframe
                    if len(comp_df) != 0: 
                        if comp_pform == 'FFld': comp_col= 'diff_FFld'
                        elif comp_pform == 'WinS': comp_col= 'diff_Wins'
                        elif comp_pform == 'LIDR': comp_col= 'diff_LiDR'
                        elif comp_pform == 'Prb1': comp_col= 'diff_Prb1'
                        elif comp_pform == 'Prb2': comp_col= 'diff_Prb2'
                        elif comp_pform == 'CoMeT1': comp_col= 'diff_CoMeT1'
                        elif comp_pform == 'CoMeT2': comp_col= 'diff_CoMeT2'
                        elif comp_pform == 'CoMeT3': comp_col= 'diff_CoMeT3'

                        comp_p_time = nearest(comp_df.index, p_time)
                        comp_lat= full_df.loc[(comp_pform, comp_p_time), 'lat']
                        comp_lon= full_df.loc[(comp_pform, comp_p_time), 'lon']
                        pform_lat = full_df.loc[(pform, p_time), 'lat']
                        pform_lon = full_df.loc[(pform, p_time), 'lon']
                        
                        # Calc dist between the 2 instruments (assuming spherical earth) outputs in km 
                        dist = mpu.haversine_distance((pform_lat, pform_lon), (comp_lat, comp_lon))
                        full_df.loc[(pform, p_time), comp_col] = dist
                    # ***
        
        ##Save to a csv file
        #This is the directory path for the output file
        full_df.to_csv(outdir+day+'_comparison_full.csv')


    #############
    if apply_dist_threshold == True:
        df = pd.read_csv(outdir+day+'_comparison_full.csv')

        col_str='pforms_within_'+str(valid_dist_diff)+'km'
        differences = df.drop(columns=['pname', 'datetime', 'lat', 'lon', 'tfast', 'Thetav', 'Thetae'])
        differences[col_str]=np.nan
        differences.drop(columns=[col_str], inplace=True)

        #figure out how many pforms are within our dist threshold for each obs (binned) time 
        differences_bool = differences.le(valid_dist_diff)
        # Add new column based on conditions and choices:
        pforms_valid= differences_bool.sum(axis=1)

        # need to subtract one since otherwise the pform will be counted as being close to itself
        df[col_str]=pforms_valid.subtract(1)
        df.to_csv(outdir+day+'_comparison_full.csv', index=False)
        ##########

    #############
    if apply_storm_threshold == True:
        for end_str in ['_full.', '.']:
            df = pd.read_csv(outdir+day+'_comparison'+end_str+'csv')
            df = df[df['datetime']<=day_cuttoff[day]]
            df.to_csv(outdir+day+'_comparison_ps'+end_str+'csv', index=False)

    #############
    if plotting == True:
        for end_str in ['_ps_full.','_full.', '_ps.', '.']:
            df = pd.read_csv(outdir+day+'_comparison'+end_str+'csv')
            df['datetime']=pd.to_datetime(df['datetime'], format='%Y-%m-%d %H:%M:%S')

            num_of_pforms=df['pname'].nunique()
            ##########
            if prox_plots == True:
                col_str='pforms_within_'+str(valid_dist_diff)+'km'
                fig, ax = plt.subplots(nrows=num_of_pforms, ncols=1, sharex=True, sharey=True, figsize=(15,10)) 
                r=0
                for label, grp in df.groupby('pname'):
                    ax[r].yaxis.set_major_locator(MultipleLocator(2))
                    ax[r].yaxis.set_minor_locator(MultipleLocator(1))
                    ax[r].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M')) #strip the date from the datetime object
                    ax[r].xaxis.set_major_locator(mdates.HourLocator())
                    ax[r].xaxis.set_minor_locator(mdates.HourLocator())
                    ax[r].grid('on', which='both', axis='both', alpha=.5)
                    ax[r].set_xlim(min(df['datetime']),max(df['datetime']))
                    grp.plot.scatter(x = 'datetime', y =col_str, color=colors[label], ax = ax[r], label = label)
                    #  grp.plot.area(x = 'datetime', y = 'pforms_within_10km', color=colors[label], ax = ax, label = label, alpha=.25, stacked=False)
                    #  grp.plot(x = 'datetime', y = 'pforms_within_10km', color=colors[label], ax = ax[r], label = label)
                    #  grp.plot(x = 'datetime', y = 'pforms_within_10km', ax = ax[r], label = label)
                    r=r+1

                plt.xlabel('Time (HH:MM UTC)')
                plt.ylabel('pforms within '+ str(valid_dist_diff))
                fig.suptitle(day+' platform proximity '+str(valid_dist_diff)+ 'km '+end_str)
                fig.tight_layout()
                plt.subplots_adjust(wspace=0, hspace=0)
                
                name= plot_config.g_TORUS_directory+'comparison/prox/'+day+'_'+str(valid_dist_diff)+'_proximity'+end_str+'png'
                print(name)
                plt.savefig(name, bbox_inches='tight')
                plt.close('all')   # close it
                fig.clf()          # no really... close it
                gc.collect()       # get tf outta here

            #########
            if comp_plots == True:
                comp_var_list =['Thetae','Thetav']
                for col_str in comp_var_list:
                    subset= True
                    if subset == True:
                        num_pform_thres = 4
                        df=df[df.loc[:,'pforms_within_'+str(valid_dist_diff)+'km']>=num_pform_thres]
                    fig,ax = plt.subplots(figsize=(15,10)) 
                    for label, grp in df.groupby('pname'):
                        ax.yaxis.set_major_locator(MultipleLocator(2))
                        ax.yaxis.set_minor_locator(MultipleLocator(1))
                        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M')) #strip the date from the datetime object
                    ax.xaxis.set_major_locator(mdates.HourLocator())
                    ax.xaxis.set_minor_locator(mdates.HourLocator())
                    ax.grid('on', which='both', axis='both', alpha=.5)
                    ax.set_xlim(min(df['datetime']),max(df['datetime']))
                    #  grp.plot.scatter(x = 'datetime', y =col_str, color=colors[label], ax = ax[r], label = label)
                    #  grp.plot.area(x = 'datetime', y = 'pforms_within_10km', color=colors[label], ax = ax, label = label, alpha=.25, stacked=False)
                    grp.plot.scatter(x = 'datetime', y = col_str , color=colors[label],linestyle='-', ax = ax, label = label)
                    #  grp.plot(x = 'datetime', y = 'pforms_within_10km', ax = ax[r], label = label)

                    plt.xlabel('Time (HH:MM UTC)')
                    plt.ylabel(col_str)
                    fig.suptitle(day+' '+col_str+' '+ end_str)
                    fig.tight_layout()
                    
                    name= plot_config.g_TORUS_directory+'comparison/comps/'+day+'_'+col_str+end_str+'png'
                    print(name)
                    plt.savefig(name, bbox_inches='tight')
                    plt.close('all')   # close it
                    fig.clf()          # no really... close it
                    gc.collect()       # get tf outta here
            ########
