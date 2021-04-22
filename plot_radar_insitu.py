#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import needed modules
######################
import matplotlib
#  matplotlib.use('agg')
#  matplotlib.use('GTK3Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.patheffects as PE
from matplotlib.collections import PatchCollection
from matplotlib.offsetbox import (AnchoredOffsetbox, DrawingArea, HPacker, TextArea)
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredAuxTransformBox
from matplotlib import ticker
from matplotlib.ticker import (LinearLocator, FixedLocator, MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.colors 
import time
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
import singledop
import pynimbus as pyn 
import tracemalloc
import shapefile
import pickle
import math 
tracemalloc.start()

import pprint
pp = pprint.PrettyPrinter(indent=4)

register_matplotlib_converters()
logging.basicConfig(filename='this_is.log')
log = logging.getLogger(__name__)
# To run with a) warnings and b) stack trace on abort
# python3 -Walways  -q -X faulthandler plot_nexrad_insitu.py

## Imports form other files
############################
if plot_config.overlays['GEO']['OX'] == True: import osmnx as ox

## Read in defns that I have stored in another file (for ease of use/consistancy accross multiple scripts)
from read_pforms import pform_names, Add_to_DATA, Platform, Torus_Insitu, Radar, Stationary_Insitu
from read_pforms import error_printing, timer, time_in_range, time_range
from radar_defns import read_from_radar_file, det_nearest_WSR, sweep_index, get_WSR_from_AWS, info_radar_files

totalcompT_start = time.time()
################################################################################################
##########
# Classes
##########
class Thermo_Plt_Vars:
    def __init__ (self, config, day, Data):
        self.Te_lab, self.Tv_lab, self.Tf_lab = "Equi. Pot. Temp [K]", "Vir. Pot. Temp [K]", "Temp (fast) [K]"
        self.Tv_GMin, self.Tv_GMax = self.find_global_max_min('Thetav', Data)
        self.Te_GMin, self.Te_GMax = self.find_global_max_min('Thetae', Data)
        self.Tf_GMin, self.Tf_GMax = self.find_global_max_min('tfast', Data)
            
        print(Data)
        if config.overlays['Colorline']['Var'] == "Thetae":
            lab, GMin, GMax = self.Te_lab, self.Te_GMin, self.Te_GMax
        elif config.overlays['Colorline']['Var'] == "Thetav":
            lab, GMin, GMax = self.Tv_lab, self.Tv_GMin, self.Tv_GMax
        elif config.overlays['Colorline']['Var'] == "tfast":
            lab, GMin, GMax = self.Tf_lab, self.Tf_GMin, self.Tf_GMax
        self.Tvar_lab, self.Tvar_GMin, self.Tvar_GMax = lab, GMin, GMax

        #settings for the thermo var being plotted on radar subplots
        ## Make a dummy plot to allow for ploting of colorbar
        self.norm =plt.Normalize(GMin, GMax)

        Z = [[0,0],[0,0]]
        levels = np.arange(GMin, GMax+1, 1)
        self.CS3 = plt.contourf(Z, levels, cmap='gnuplot')
        #  self.CS3 = plt.contourf(Z, levels, cmap=cmocean.cm.thermal)
        plt.clf()

    # * * *
    def find_global_max_min(self, var, Dict):
        '''determine the global max and min across all platforms for a given var
        '''
        val_hold = []
        for p in Dict.values():
            if var == 'Thetav':
                if hasattr(p, 'Tv_Min') == True:  val_hold.append(p.Tv_Min)
                if hasattr(p, 'Tv_Max') == True:  val_hold.append(p.Tv_Max)
            if var == 'Thetae':
                if hasattr(p, 'Te_Min') == True:  val_hold.append(p.Te_Min)
                if hasattr(p, 'Te_Max') == True:  val_hold.append(p.Te_Max)
            if var == 'tfast':
                if hasattr(p, 'Tf_Min') == True:  val_hold.append(p.Tf_Min)
                if hasattr(p, 'Tf_Max') == True:  val_hold.append(p.Tf_Max)
        if len(val_hold) != 0: global_min, global_max = min(val_hold), max(val_hold)
        return global_min, global_max

#### * * * * * * * * * * * * * * * * * * * * * * * * * *
### set up Master_Plt class (along with subclasses R_Plt (radar plots), and TS_Plt (timeseries plot))
class Master_Plt:
    def __init__(self, Data, config):
        self.Data = Data #(the data dict)
        self.config = config

        ## Set some parameters for the plot such as fontsize etc
        #### * * * * * * * * * * * * * * * * * * * * * * * * * *
        #  print(plt.rcParams)
        plt.rc('font', size= 15)      # controls default text sizes
        plt.rc('axes', labelsize= 21) # fontsize of the axes title, and x and y labels
        plt.rc('legend', fontsize= 23, borderpad=.5, facecolor='white', 
               edgecolor= 'black', shadow=True, fancybox=True, framealpha=1)# legend fontsize
        plt.rc('figure', titlesize= 50, facecolor='white')  # fontsize of the figure title
        self.leg_title_font={'size':25, 'weight':'bold'}
        self.Radar_title_font= {'fontsize':15}
        #  self.Radar_title_font= {'fontsize':40}

        ## Establish the Plot size and layout
        #### * * * * * * * * * * * * * * * * *
        def radar_layout(config):
            #if you specify a start or end time it will be assigned here otherwise will be set to none (full data set)
            if config.radar_controls['layout']['Rows'] == None:
                R_rows = math.trunc(len(config.r_mom) / 2)
            else:
                R_rows = config.radar_controls['layout']['Rows']

            if config.radar_controls['layout']['Cols'] == None:
                if ((len(config.r_mom) / 2) % 2) in [0,1]: R_cols = 2
                else: R_cols = 3
            else: R_cols = config.radar_controls['layout']['Cols']
            
            return R_rows, R_cols

        self.title_spacer=.93
        # This is the layout for radar and time series subplots
        if len(config.Time_Series) != 0 and len(config.r_mom) != 0:
            self.fig= plt.figure(figsize=(32,20))
            if len(config.Time_Series) == 1:
                outer_hratio = [4,2] #the height ratio between the space alocated for radar imagaes and timeseries
                ts_hratio=[1] #the relative height alocated to each timeseries [within the space allocated for all timeseries]
            if len(config.Time_Series) == 2:
                outer_hratio = [4,3] 
                if 'Wind' in config.Time_Series: ts_hratio=[2,3]
                else: ts_hratio=[1,1]
            self.outer_gs= GridSpec(nrows=2, ncols=1, height_ratios=outer_hratio, hspace=.1)

            self.ts_gs = GridSpecFromSubplotSpec(len(config.Time_Series), 1, subplot_spec=self.outer_gs[1, :], height_ratios=ts_hratio, hspace=0)
            self.R_rows, self.R_cols = radar_layout(config)
            self.r_gs = GridSpecFromSubplotSpec(self.R_rows, self.R_cols, subplot_spec=self.outer_gs[0, :])

        # This is the layout for time series only
        if len(config.Time_Series) != 0 and len(config.r_mom) == 0:
            print('this is the layout for time series only')
            #  self.fig= plt.figure(figsize=(32,5))
            self.fig= plt.figure(figsize=(32,10))
            if len(config.Time_Series) == 1:  ts_hratio=[1]
            if len(config.Time_Series) == 2:
                if 'Wind' in config.Time_Series: ts_hratio=[2,3]
                else: ts_hratio=[1,1]
            self.ts_gs = GridSpec(len(config.Time_Series), 1, height_ratios=ts_hratio, hspace=0)

        # This is the layout for radar only
        if len(config.r_mom) != 0 and len(config.Time_Series) == 0:
            plt.rc('font', size= 10)      # controls default text sizes
            plt.rc('figure', titlesize= 25, facecolor='white')  # fontsize of the figure title
            plt.rc('axes', labelsize= 10) # fontsize of the axes title, and x and y labels
            print('this is the layout for radar plots only')
            #  self.fig= plt.figure(figsize=(15,25))

            self.R_rows, self.R_cols = radar_layout(config)
            if self.R_cols == 3: 
                self.fig= plt.figure(figsize=(30,10))
            if self.R_cols == 2: 
                self.fig= plt.figure(figsize=(20,20.5))
            if self.R_rows == 1: 
                self.title_spacer=.68
            if self.R_rows == 2: 
                self.title_spacer=.8
            self.Radar_title_font= {'fontsize':20}
            self.r_gs = GridSpec(nrows=self.R_rows, ncols=self.R_cols, hspace=0)

        ## Establish some vars that will be helpful to the plot later on
        #### * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        #if we are plotting radar
        if len(config.r_mom) != 0:
            # The extent of the area to be plotted
            self.Domain = Platform.getLocation(Data[config.radar_controls['Centered_Pform']], scan_time=Data['P_Radar'].Scan_time, 
                                               offsetkm= config.radar_controls['offsetkm'][config.Radar_Plot_Type])
            #  self.Domain_Bbox = Bbox.from_extents(self.Domain.xmin, self.Domain.ymin, self.Domain.xmax, self.Domain.ymax)

            # Define pyart display object for plotting radarfile
            self.display = pyart.graph.RadarMapDisplay(Data['P_Radar'].rfile)

            # Set the projection of the radar plot, and transformations
            self.R_Proj = self.display.grid_projection
            self.Proj = ccrs.PlateCarree()

    # * * * * * * *
    def tick_grid_settings(self, ax, scan_time, radar=None, ts=None, interval=None, twinax=False):
        ##Axis Limits
        # # # # # # # #
        #if you are calling this for a timeseries subplot set the limits for the plot
        if ts != None:
            # Xaxis
            # if desired this subsets the timerange that is displayed in the timeseries
            if scan_time != None and self.config.tseries_control['Axis_Extents']['X'] != None:
                ax.set_xlim(scan_time - timedelta(minutes=self.config.tseries_control['Axis_Extents']['X']),
                            scan_time + timedelta(minutes=self.config.tseries_control['Axis_Extents']['X']))

            # Yaxis
            if ts in ['Thetav','Thetae']:
                if self.config.tseries_control['Axis_Extents']['Y'] != None: ax.set_ylim(335,350)
                else: ax.set_ylim(TVARS.Tvar_GMin, TVARS.Tvar_GMax)

            #if plotting the wind time series and its the first ax being det
            elif ts == 'Wind' and twinax == False: ax.set_ylim(0)
            #if plotting the wind time series and its the twin ax being det
            elif ts == 'Wind' and twinax == True: ax.set_ylim(0, 360)

        # + + + + + + + + + + + + ++ + +
        ##Tick Locator
        # # # # # # # #
        #if you are calling this for a timeseries subplot locate the ticks
        if ts != None:
            # Xaxis
            ax.xaxis.set_major_locator(mdates.AutoDateLocator())
            # set up minor ticks (should be a multiple of ten intervals (ie 10,20,30... min spans)
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))

            # Yaxis
            if ts in ['Thetav','Thetae']:
                # set up major tick marks (this is set up to go by 5's will want to change for diff vars)
                ax.yaxis.set_major_locator(MultipleLocator(5))
                # set up minor ticks (this have it increment by 1's will want to change for diff vars)
                ax.yaxis.set_minor_locator(AutoMinorLocator(5))
            elif ts == 'Wind' and twinax == False:  ax.yaxis.set_major_locator(LinearLocator(numticks=5))
            elif ts == 'Wind' and twinax == True:  ax.yaxis.set_major_locator(FixedLocator(np.arange(0, 450, 90)))

        #if you are calling this for a radar subplot locate the ticks
        if radar != None and interval != None:
            def calc_ticks(a,b,interval):
                a_rounded, b_rounded = round(a/interval)*interval, round(b/interval)*interval
                ticks = np.arange(a_rounded-interval, b_rounded+interval, interval)
                return ticks
            #Det the boundaries of the domain
            x0, x1 = ax.get_xbound()
            y0, y1 = ax.get_ybound()

            #Det the tick locations
            xticks, yticks =  calc_ticks(x0, x1, interval), calc_ticks(y0, y1, interval)
            #set the calculated ticks on the graph
            ax.set_xticks(xticks)
            ax.set_yticks(yticks)
            ax.locator_params(nbins=6)

            #  ax.xaxis.set_major_locator(MaxNLocator(6, steps=[1, 5, 10]))
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            #  ax.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 5, 10]))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))

            #set the calculated ticks on the graph
            plt.setp(ax.get_xticklabels(), visible=True)
            plt.setp(ax.get_yticklabels(), visible=True)

        # + + + + + + + + + + + + ++ + +
        ##Tick Display Characteristics
        # # # # # # # # # # # # # # # #
        ax.tick_params(which='minor', axis='both', width=2, length=7, color='grey')
        ax.tick_params(which='major', axis='both', width=2, length=10, color='black')
        if radar != None:
            ax.tick_params(which='both', axis='both', direction='in', top=True, bottom=False, labeltop=True, labelbottom=False)
            ax.tick_params(which='both', axis='y', labelrotation=90)

        # + + + + + + + + + + + + ++ + +
        ##Grid Display Characteristics
        # # # # # # # # # # # # # # # #
        if ts != None:
            ax.tick_params(which='major', axis='both', grid_linewidth=2.5)
            ax.tick_params(which='major', axis='y', grid_color='grey', grid_linewidth=2, grid_linestyle='--', grid_alpha=.8)

        if radar != None:
            ax.tick_params(which='both', axis='both', grid_linewidth=2, grid_color='grey', grid_alpha=.3)
        ax.tick_params(which='minor', axis='both', grid_linestyle=':', grid_alpha=.3)
        ax.grid(which='both')
        ax.set_axisbelow('line')

        # + + + + + + + + + + + + ++ + +
        ##Tick Label Formatter
        # # # # # # # #
        if ts != None:
            # Xaxis
            ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M')) #strip the date from the datetime object
            # Yaxis
            if ts == 'Wind' and twinax == False: ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))

        if radar != None:
            # Add km to x-axis.
            def km(x, pos):
                x = int(x/1000)
                return ' {} km'.format(x)
            #since this is square grid both the x and y formatter will be the same
            ax.yaxis.set_major_formatter(FuncFormatter(km))
            ax.xaxis.set_major_formatter(FuncFormatter(km))
        # + + + + + + + + + + + + ++ + +

    # * * * * * * *
    def time_series(self, ts, ax_n, fig, TVARS, Data):
        ''' Plot a time series (sub)plot; could be of pvar info from various intstruments or of wind info
        ----
        INPUTS: ts: ........fill in
                Data: dict as described in plotting defn
                print_long & e_test: bool str as described in the ppi defn '''
        if self.config.print_long == True: print('~~~~~~~~~~~Made it into time_series~~~~~~~~~~~~~~~~~')
        #  t_TS.T_Plt_settings(ts,ax_n)
        TSleg_elements = [] #empty list to append the legend entries to for each subplot that is actually plotted
        ## MAKE THE TIMESERIES
        #### * * * * * * * * *
        if ts == 'histogram':
            hist, bins = np.histogram(Data['P_Radar'].rfile.fields[self.config.tseries_control['H']['var']]['data'][~np.isnan(Data[
                'P_Radar'].rfile.fields[self.config.tseries_control['H']['var']]['data'])], bins=150)
            bins = (bins[1:]+bins[:-1])/2.0
            ax_n.plot(bins, hist)
            ax_n.set_xlim(self.config.tseries_control['H']['xmin'], self.config.tseries_control['H']['xmax'])
            #  ax_n.set_yscale('log')
            #  ax_n.axvline(x=3,c='k')
            ax_n.set_xlabel(self.config.tseries_control['H']['var'])
            ax_n.set_ylabel('Count')

        elif ts in ['Thetae', 'Thetav', 'Wind']:
            if ts in ['Thetav', 'Thetae']:
                for p in Data.values():
                    if isinstance(p, Torus_Insitu):
                        if self.config.print_long == True: print('Plotting '+str(p.name)+' on time series')
                        ## If the platform matches a type listed for masking it uses masked data for the time series; 
                        #  else use the unmasked data
                        if p.type in self.config.tseries_control['Mask'][:]:
                            if ts == "Thetae":   plotting_data= p.Thetae_ts_mask_df
                            elif ts == "Thetav":   plotting_data= p.Thetav_ts_mask_df
                        else: plotting_data= p.df[ts].values

                        ## Plot
                        #assigning label= is what allows the legend to work
                        ax_n.plot(p.df['datetime'], plotting_data, linewidth=3, color=p.l_color, label=p.leg_str)

                        TSleg_entry = Line2D([], [], label=p.leg_str, linewidth=12, color=p.l_color)
                        TSleg_elements.append(TSleg_entry)
                ## Set up XY axes tick locations and Labels
                if len(self.config.r_mom) != 0: scan= Data['P_Radar'].Scan_time
                else: scan= None 
                self.tick_grid_settings(ax=ax_n, ts=ts, scan_time=scan)

                if ts == "Thetae":   ylab = TVARS.Te_lab
                elif ts == "Thetav":   ylab = TVARS.Tv_lab
                ax_n.set_ylabel(ylab)


            # * * *
            if ts == 'Wind':
                p = Data[self.config.tseries_control['Wind_Pform']]
                if self.config.print_long == True: print('Plotting '+str(p.name)+' on time series')

                ax_n.plot(p.df['datetime'], p.df['spd'], label= 'Wind Spd')
                ax_n.fill_between(p.df['datetime'], p.df['spd'], 0)
                ax_n.set_ylabel('Wind Speed (kn)')
                self.tick_grid_settings(ax=ax_n, ts=ts,scan_time=Data['P_Radar'].Scan_time)
                TSleg_entry = Line2D([], [], label='Wind Spd', linewidth=12, color='tab:blue')
                TSleg_elements.append(TSleg_entry)

                ax_2 = ax_n.twinx()
                ax_2.plot(p.df['datetime'], p.df['dir'], '.k', linewidth=.05, label='Wind Dir')
                ax_2.set_ylabel('Wind Dir ($^{\circ}$)')
                self.tick_grid_settings(ax=ax_2, ts=ts, twinax=True, scan_time=Data['P_Radar'].Scan_time)
                TSleg_entry = Line2D([], [], marker='.', color='black', label='Wind Dir', markersize=26)
                TSleg_elements.append(TSleg_entry)

            ## Plot legend
            leg = ax_n.legend(handles= TSleg_elements, scatterpoints=3, loc='center left')
            if ts == 'Wind':
                leg.set_title(self.config.tseries_control['Wind_Pform'], prop=self.leg_title_font)
                leg.remove()
                ax_2.add_artist(leg)

            #if plotting more than one timeseries only include the xaxis label & ticks to the bottom timeseries
            ax_n.label_outer()

            #if you are not making the bottom timeseries
            if ax_n.get_subplotspec().rowspan.stop != len(self.config.Time_Series):
                ax_n.spines['bottom'].set_linewidth(5)
                # if you are doing two thermo plots while also plotting radar remove one of the legends (for space)
                if (len(self.config.r_mom) != 0) & (ts in ['Thetav', 'Thetae']): leg.remove()

            ## If makeing the timeseries in conjunction with radar subplots set up vertical lines that indicate the time of
            #  radar scan and the timerange ploted (via filled colorline) on the radar plots
            if len(self.config.r_mom) != 0:
                ax_n.axvline(Data['P_Radar'].Scan_time, color='r', linewidth=4, alpha=.5)
                ax_n.axvspan(Data['P_Radar'].Scan_time - timedelta(minutes=self.config.overlays['Colorline']['cline_extent']),
                             Data['P_Radar'].Scan_time + timedelta(minutes=self.config.overlays['Colorline']['cline_extent']), 
                             facecolor='0.5', alpha=0.4)
            ##Label Formatter
            # # # # # # # # #
            ax_n.set_xlabel('Time (UTC)')
            ax_n.yaxis.set_label_coords(-.03, .5)
            #  make_axes_area_auto_adjustable(ax_n)

            ###################
        if self.config.print_long == True: print('~~~~~~~~~~~Made it through time_series~~~~~~~~~~~~~~')
    
    # * * * * * *  *
    def radar_subplots(self, mom, day, ax_n, fig, TVARS, Data, leg):
        ''' Plots each of the radar subplots including the marker overlays corresponding to the additional platforms
        ----
        INPUTS: mom: str, indicates which radar moment you would like to plot (ie Reflectivity, Velocity, Spectrum Width etc)
                Data: dict as described in the plotting defn
                leg: bool str whether or not you want a legend associated with this particular subplot
        '''
        if self.config.print_long == True: print('~~~~~~~~~~~made it into radar_subplots~~~~~~~~~~~~~~')


        ## SET UP VARS FOR EACH RADAR MOMENTS
        #####################################
        testing_testing= False
        testing_bound=None
        if mom in ['refl','refl_unmasked', 'refl_despeck']:
            c_scale = 'pyart_HomeyerRainbow'
            if Data['P_Radar'].name == 'WSR88D':
                sweep = Data['P_Radar'].swp[0]
                vminb, vmaxb = -10., 75.
            else:
                sweep = Data['P_Radar'].swp
                if Data['P_Radar'].name in pform_names('KA'): vminb, vmaxb= -30., 30.
                elif Data['P_Radar'].name == 'NOXP': vminb, vmaxb= -10, 60
            #****
            if mom == 'refl':
                c_label= 'Reflectivity [dbz]' 
                if Data['P_Radar'].name == 'WSR88D': field = 'reflectivity'
                else: field= 'refl_fix'
            if mom == 'refl_unmasked':
                c_label= 'Reflectivity [dbz]' 
                if Data['P_Radar'].name == 'WSR88D': field = 'reflectivity'
                elif Data['P_Radar'].name in pform_names('KA'): field= 'reflectivity'
                elif Data['P_Radar'].name =='NOXP': field= 'DBZ'
            elif mom == 'refl_despeck':
                c_label= 'Reflectivity DS [dbz]'
                field= 'refl_despeck'

        ###
        elif mom in ['vort', 'vort_smooth']:
            if mom == 'vort_smooth':
                c_label = 'Vort smooth 'r' [$\frac{1}{s}$]'
                field = 'vort_smooth'
            elif mom == 'vort':
                c_label = 'Vort 'r' [$\frac{1}{s}$]'
                field = 'vort'
            c_scale =  cmocean.cm.curl
            vminb,vmaxb, sweep = -.15, .15,  Data['P_Radar'].swp
        ###
        elif mom in ['vel', 'vel_unmasked', 'vel_unfixed', 'vel_despeck', 'vel_smooth', 'vel_savgol']: 
            c_scale = cmocean.cm.balance
            if Data['P_Radar'].name == 'WSR88D':
               vminb, vmaxb, sweep = -40., 40., Data['P_Radar'].swp[1]
            else:
               vminb, vmaxb, sweep = -30., 30., Data['P_Radar'].swp
               #  vminb, vmaxb, sweep = -30., 30., Data['P_Radar'].swp
            #****
            if mom == 'vel':
                c_label = 'Radial Velocity [m/s]'
                if Data['P_Radar'].name == 'WSR88D': field= 'velocity'
                else: field= 'vel_fix'
            elif mom == 'vel_unmasked':
                c_label = 'Radial Velocity [m/s]'
                if Data['P_Radar'].name == 'WSR88D': field= 'velocity'
                elif Data['P_Radar'].name in pform_names('KA'): field= 'corrected_velocity'
                elif Data['P_Radar'].name =='NOXP': field= 'corrected_velocity'
            elif mom == 'vel_unfixed':
                c_label= 'Velocity Unfixed [m/s]'
                if Data['P_Radar'].name in pform_names('KA'): field= 'velocity'
                elif Data['P_Radar'].name == 'NOXP': field = 'VEL' 
            elif mom == 'vel_despeck':
                c_label= 'Radial Velocity DS [m/s]'
                field= 'vel_despeck'
            elif mom == 'vel_smooth':
                c_label= 'Radial Velocity Smoothed [m/s]'
                field= 'vel_smooth'
            elif mom == 'vel_savgol':
                c_label= 'Radial Velocity Savgol [m/s]'
                field= 'vel_savgol'

        ###
        elif mom == 'specw':
            c_label, c_scale = 'Spectrum Width 'r' [$\frac{m}{s}$]', 'cubehelix_r'
            field, vminb, vmaxb, sweep = 'specw_fix', 0., 4, Data['P_Radar'].swp

        ###
        elif mom in ['vel_texture','vel_texture_smoothed', 'vel_texture_dea', 'vel_test']:
            c_scale = 'bone_r'#'PuBuGn'#'BuPu'#cmocean.cm.dense#cmocean.cm.speed#Greys #'pyart_balance'
            #  c_scale = 'Greys'#'pyart_balance'
            vminb, vmaxb, sweep =  0., 3, Data['P_Radar'].swp
            #  vminb, vmaxb, sweep =  0., 10, Data['P_Radar'].swp
            if mom == 'vel_texture':
                c_label = 'Radial Velocity Texture (no DEA) STD Dev' 
                field = 'vel_texture'
            if mom == 'vel_texture_smoothed':
                c_label = 'Radial Velocity Texture Smoothed STD Dev' 
                field = 'vel_texture_smoothed'

            elif mom == 'vel_texture_dea':
                c_label = 'Radial Velocity Texture (DEA) STD Dev '
                field= 'vel_texture_dea'
            elif mom == 'vel_test': 
                #  testing_testing= False
                #  testing_bound=None
                testing_testing= True
                testing_bound=[0,.5,1,1.5,2,2.5,3,3.5]
                c_label= 'Radial Velocity masking test'
                field= 'vel_test'

        ###
        elif mom in ['vel_grad', 'vel_savgol_axis1', 'vel_savgol_axis0']:
            sweep = Data['P_Radar'].swp
            c_scale = cmocean.cm.balance  
            if mom == 'vel_grad':
                p_title, c_label= 'Vel Gradient', 'Gradient m/s'
                field, vminb, vmaxb= 'vel_gradient', -5., 5
            elif mom == 'vel_savgol_axis1':
                c_label= 'Vel Savgol axis1 '
                field, vminb, vmaxb= 'vel_savgol_axis1', -1., 1
            elif mom == 'vel_savgol_axis0':
                c_label= 'Vel savgol axis0'
                field, vminb, vmaxb= 'vel_savgol_axis0', -2., 2
            elif mom == 'vel_savgol_grad':
                p_title, c_label= 'Radial Gradient (savgol smoothed)', 'm/s'
                field, vminb, vmaxb= 'vel_savgol_grad', -5., 5
            elif mom == 'vel_grad0_0':
                p_title, c_label= 'Radial Velocity (savgol 0) A0', 'Velocity [m/s]'
                field, vminb, vmaxb= 'vel_gradient0_0', -30., 30
            elif mom == 'vel_grad1_0':
                p_title, c_label= 'Radial Acceleration (savgol 1) A0', 'm/s^2'
                field, vminb, vmaxb= 'vel_gradient1_0', -2., 2
            elif mom == 'vel_grad2_0':
                p_title, c_label= 'Radial 3rd Deriv (savgol 2) A0', 'm/s^3'
                field, vminb, vmaxb= 'vel_gradient2_0', -0.2, 0.2
            elif mom == 'vel_grad0_1':
                p_title, c_label= 'Radial Velocity (savgol 0) A1', 'Velocity [m/s]'
                field, vminb, vmaxb= 'vel_gradient0_1', -30., 30
            elif mom == 'vel_grad1_1':
                p_title, c_label= 'Radial Acceleration (savgol 1) A1', 'm/s^2'
                field, vminb, vmaxb= 'vel_gradient1_1', -1., 1
            elif mom == 'vel_grad2_1':
                p_title, c_label= 'Radial 3rd Deriv (savgol 2) A1', 'm/s^3'
                field, vminb, vmaxb= 'vel_gradient2_1', -0.1, 0.1
            elif mom == 'sim_vel': 
                p_title, c_label, c_scale = 'Sim', 'Difference', 'Greys'#'pyart_balance'
                field, vminb, vmaxb= 'difference', 0., 2.

        ###
        elif mom in ['az_shear', 'div_shear']:
            sweep = Data['P_Radar'].swp
            if mom == 'az_shear':
                c_label, c_scale = 'Azimulathal Shear (s$^{-1}$)', 'RdBu'
                field, vminb, vmaxb = 'az_shear', -15, 15 
        ###
        elif mom == 'diff':
            p_title, c_label, c_scale = 'Difference (fixed - original)', 'Difference', cmocean.cm.balance#'pyart_balance'
            field, sweep = 'difference', Data['P_Radar'].swp

            diff_field = Data['P_Radar'].rfile.get_field(Data['P_Radar'].swp, 'difference', copy=False)
            #  diff_field = Data['P_Radar'].rfile.fields['difference']['data']
            vminb, vmaxb = np.min(diff_field), np.max(diff_field)
        else: print("Hey what just happened!\n Check the Radar moments for spelling")
        #####################################
        #  print("This is the id for "+ field+' :'+str(id(Data['P_Radar'].rfile.get_field(Data['P_Radar'].swp, field, copy=False))))

        ## set colors to highlight where the max/min are exceded
        colormap = copy.copy(matplotlib.cm.get_cmap(name=c_scale))
        colormap.set_under('aqua')
        colormap.set_over('magenta')

        ## Plot the radar
        #  ax_n.set_title(p_title, y=-.067, fontdict=self.Radar_title_font)
        self.display.plot_ppi_map(field, sweep, ax=ax_n, cmap=colormap, vmin=vminb, vmax=vmaxb,
                                  width=self.config.radar_controls['offsetkm'][self.config.Radar_Plot_Type]*2000,
                                  height=self.config.radar_controls['offsetkm'][self.config.Radar_Plot_Type]*2000, 
                                  title_flag=False, colorbar_flag=False, embelish=False)
        
        # Has to be here or it doesn't work
        ax_n.set_extent(self.Domain)
        self.tick_grid_settings(ax=ax_n, radar=True, scan_time=None, interval=.5*1000)
        #  self.tick_grid_settings(ax=ax_n, radar=True, scan_time=None, interval=10*1000)
        ax_n.grid(True)

        # Do it again here to restore extent
        ax_n.set_extent(self.Domain)
        

        # * * *
        ## PLOT PLATFORMS AS OVERLAYS(ie marker,colorline etc) ON RADAR
        def plot_markers(p, lon, lat, txt_label, lab_shift=(.009,.002), lab_c='xkcd:pale blue'):
            '''Plot a marker at the locations of the various platforms: Optional text labels on plot as well'''
            if isinstance(p, Stationary_Insitu):
                i=0
                for x, y in zip(lon, lat):
                    ax_n.plot(x, y, transform=self.Proj, marker=p.m_style, mfc=p.m_color, linestyle=None,
                              ms=p.m_size, mec='k', label=p.leg_str if i==0 else "", zorder=10) 
                    i=i+1
            else:    
                ax_n.plot(lon, lat, transform=self.Proj, marker=p.m_style, mfc=p.m_color, linestyle=None,
                          ms=p.m_size, mec='k', label=p.leg_str, zorder=10, path_effects=[PE.withStroke(linewidth=4, foreground='k')])

            ## Optional textlabel on plot
            if self.config.overlays[p.type]['Label'] == True:
                if isinstance(txt_label, str):
                    ax_n.text(lon+lab_shift[0], lat-lab_shift[1], txt_label, transform=self.Proj,
                              path_effects=[PE.withStroke(linewidth=4, foreground=lab_c)])
                if isinstance(txt_label, pd.core.series.Series):
                    for x, y, lab in zip(lon, lat, txt_label):
                        ax_n.text(x+lab_shift[0], y-lab_shift[1], lab, transform=self.Proj,
                                  path_effects=[PE.withStroke(linewidth=4, foreground=lab_c)])

        def plot_TORUSpform(self, Data, p, border_c='xkcd:light grey'):
            '''Plot the filled line that indicates the thermo values on the Rplt'''
            #grab the subset of data of +- interval around radar scan
            p_sub, p_deploy = p.grab_pform_subset(p, Data, time_offset=self.config.overlays['Colorline']['cline_extent'])

            #if there is data for the platform that falls within the time and location of interest
            if p_deploy == True:
                ##READ in data from p_sub
                x, y, C = p_sub['lon'].values, p_sub['lat'].values, p_sub[self.config.overlays['Colorline']['Var']].values

                ##PLOT the line that extends +/- min from the platform location;
                if self.config.overlays['Colorline'][p.type] == True:
                    #  The colorfill indicates values of the specifiec p_var (ie Thetae etc)
                    ax_n.scatter(x, y, c=C, cmap=cmocean.cm.thermal, norm=TVARS.norm, alpha=.7, marker='o',
                                 s=7.5, zorder=3, transform=self.Proj)
                    #  background of the colorline
                    ax_n.plot(x, y, c=border_c, lw=9, transform=self.Proj)
                    #  plot a dot at the end of the colorline in the direction the platform is moving
                    ax_n.plot(x[-1], y[-1], color='k', marker='.', markersize=9, zorder=4, transform=self.Proj)

                    ##PLOT windbarbs
                    #  p_sub.iloc[::x,col_index] returns every x'th value
                    stationplot = StationPlot(ax_n, p_sub.loc[::30, 'lon'], p_sub.loc[::30, 'lat'], transform=self.Proj)
                    stationplot.plot_barb(p_sub.loc[::30, 'U'], p_sub.loc[::30, 'V'], sizes=dict(emptybarb=0), length=7)

                ##DET the row of the pform data at scantime (for plotting marker)
                C_Point=p_sub.loc[p_sub['datetime'].sub(Data['P_Radar'].Scan_time).abs().idxmin()]
                return C_Point.lon, C_Point.lat, True

            elif p_deploy == False:
                if self.config.print_long == True: print('The platform was not deployed at this time')
                return False, False, False

        def rhi_spokes_rings(pform, ax= ax_n):
            ''' Plot the RHI spoke and ring for a radar '''
            # if there is not actually rhi info then it will not plot a ring and not stop the code
            if np.isnan(pform.rhib) == True or np.isnan(pform.rhie)== True: print('Could not plot RHI spokes')
            #produce spoke and ring
            else:
                latArray, lonArray = [], []
                if pform.type == 'KA': radius= 20.905
                for bearing in range(int(pform.head + pform.rhib), int(pform.head + pform.rhie+1), 10): #degrees of sector
                    lat2, lon2 = Platform.getLocation(pform, radius, scan_time=Data['P_Radar'].Scan_time, given_bearing=bearing)
                    #plot spokes
                    ax.plot([pform.lon, lon2], [pform.lat, lat2], color='k', linewidth=0.5, linestyle=":", transform=self.Proj)
                    #append to array for connecting circle
                    latArray.append(lat2)
                    lonArray.append(lon2)
                #plot the circle that connects the spokes
                ax.plot(lonArray, latArray, c='k', lw=.5, transform=self.Proj)

        # * * *
        #  iterate over each object contained in dict Data (returns the actual objects not their keys)
        for p in Data.values():
            if self.config.print_long == True: print(p.name)
            ######
            #Plot Inistu Torus Platforms (if desired and available)
            if isinstance(p, Torus_Insitu):
                LON, LAT, valid_sites = plot_TORUSpform(self, Data, p, border_c='k')
                txt_label = p.name
                
            #Plot Stationary Inistu Platforms (aka mesonets and ASOS) if desired and available
            elif isinstance(p, Stationary_Insitu):
                # determine if any of the sites fall within the plotting domain
                sites_subdf, valid_sites = p.grab_pform_subset(p, Data, bounding= self.Domain)
                LON, LAT = sites_subdf.lon, sites_subdf.lat
                if len(sites_subdf) !=0:
                    if p.name == 'OK_M': txt_label = sites_subdf.stid
                    else: txt_label = sites_subdf.Stn_ID
            
            #Plot Radar Platforms (if desired and available)
            elif isinstance(p, Radar):
                if p.type == 'MAINR':
                    valid_sites=False
                elif p.type != 'MAINR':
                    #det if radar(s) are loc within the plot domain at the time of the plot
                    sites_subdf, valid_sites = p.grab_pform_subset(p, Data, bounding= self.Domain)
                    if p.type == 'WSR':
                        LON, LAT, txt_label = sites_subdf.lon, sites_subdf.lat, sites_subdf
                    #if these conditions are met then plot the radar marker(s)
                    if p.type == 'KA' or p.type == 'NOXP':
                        LON, LAT, txt_label = p.lon, p.lat, p.name
                        ## Plot RHI spokes
                        if p.type == 'KA' and self.config.overlays['KA']['RHI_ring'] == True: 
                            rhi_spokes_rings(p)

            ##PLOT the platform marker
            ##########################
            # if there are sites within the domain plot the markers and include in the legend
            if self.config.overlays[p.type]['Marker'] == True and valid_sites == True: 
                plot_markers(p, LON, LAT, txt_label)

        # * * *
        ## PLot LSR and Tor damage paths
        if self.config.overlays['LSR']['Marker']== True:
            link = "https://www.spc.noaa.gov/climo/reports/19"+day[4:6]+day[6:8]+"_rpts_filtered.csv"
            hazard_colors={'tornado': 'red', 'hail':'green', 'wind': 'blue'}
            for hazard in ['tornado', 'hail', 'wind']:
                lsr_reports = pyn.get_spc_storm_reports(link, type_of_df= hazard)
                lat, lon =lsr_reports.points.lat_lon_to_plot()
                ax_n.scatter(lon, lat, color=hazard_colors[hazard], transform=self.Proj)
        
        if self.config.overlays['Tracks']['Marker']== True:
            # to Download Data -> https://apps.dat.noaa.gov/StormDamage/DamageViewer/
            # enter dates, and zoom into wanted tracks
            # click "extract toolbox" and follow instructions # click shapefile
            # download the created shapefile file

            # once you download data, you should have 4 files with the same name but with endings .dbf, .prj, .shp, .shx
            # include the name, but not any file ending in the directory path
            #  track_data = self.config.g_TORUS_directory+day+'/data/reports/extractDamage/nws_dat_damage_paths'
            #  track_data = self.config.g_TORUS_directory+day+'/data/reports/extractDamage/nws_dat_damage_polys'
            track_data = self.config.g_TORUS_directory+day+'/data/reports/extractDamage/nws_dat_damage_pnts'

            # pull the data out into ShapeRecs object
            data  = shapefile.Reader(track_data)
            shapeRecs = data.shapeRecords()

            # to look at file, shapeRecs[x].record (x is any number you want) will print out info (date, tor rating, etc.)
            # each shapeRecs record represents a single tornado track
            # to get the actual data, shapeRecs[x].shape.points will give you an array filled with lon/lat tuples

            # create a DataFrame with tornado date, rating, and track locations
            # run shapeRecs[0].record to see if there is any other data you want to include in the DataFrame
            tor_tracks = pd.DataFrame()
            track_info= {}
            for track in np.arange(0, len(shapeRecs)):
                track_info['date'] = shapeRecs[track].record[0]
                track_info['rating'] = shapeRecs[track].record[11]
                track_info['lons'] = [shapeRecs[track].shape.points[i][0] for i in np.arange(0,len(shapeRecs[track].shape.points))]
                track_info['lats'] = [shapeRecs[track].shape.points[i][1] for i in np.arange(0,len(shapeRecs[track].shape.points))]
                tor_tracks = tor_tracks.append(track_info, ignore_index=True)
            #  print(shapeRecs[0].record)
            #  print(shapeRecs)
            #  print(len(shapeRecs))  #  print(shapeRecs.record)
            #  print('999999999')
            # quick example plot to make sure data pulled correctly 
            for tor in np.arange(0, len(tor_tracks)):
                ax_n.scatter(tor_tracks.loc[tor]['lons'], tor_tracks.loc[tor]['lats'], label=tor_tracks.loc[tor]['rating'], transform=self.Proj)
                #  ax_n.plot(tor_tracks.loc[tor]['lons'], tor_tracks.loc[tor]['lats'], label=tor_tracks.loc[tor]['rating'], transform=self.Proj)
                #  ax_n.plot(tor_tracks.loc[tor]['lons'], tor_tracks.loc[tor]['lats'], transform=self.Proj)
            #  ax_n.legend()
            #  print(shapeRecs[4].record[2])
            #shapeRecs[4]['lats']
            #  tortrack = pickle.loadopen('20190314.p', 'rb'))
            #  tor_num=4
            #  lats, lons = tortrack.loc[tor_num]['lats'],tortrack.loc[tor_num]['lons']
            
            #  for tor_event in tor_dict:
                #  tortrack = pickle.load(open('{}.p'.format(tor_event), 'rb'))
                #  for SN, tor_num in zip(tor_dict[tor_event][0],tor_dict[tor_event][1]):
                    #  lats, lons = tortrack.loc[tor_num]['lats'],tortrack.loc[tor_num]['lons']
                    #  snlat, snlon = probe_locs[SN][0],probe_locs[SN][1]
                    #  xs, ys = [],[]
                    #  for i,val in enumerate(lats):
                        #  temp = get_dist(lons[i], lats[i],snlon, snlat)
                        #  xs.append(temp[0])
                        #  ys.append(temp[1])
                    #  plt.plot(np.asarray(xs)/1000, np.asarray(ys)/1000, color=colors[icol], label=tor_event)
                #  icol+=spyi1
        #  print('HHHHHHHHHHHHHHHHH')
        if self.config.overlays['Contour']['Lines']==True:
            unsmoothed_contourdata = Data['P_Radar'].rfile.get_field(sweep, self.config.overlays['Contour']['Var'])
            # smooth out the lines
            contourdata = ndimage.gaussian_filter(unsmoothed_contourdata, sigma=2)

            lat, lon, _ = Data['P_Radar'].rfile.get_gate_lat_lon_alt(sweep)
            levels = np.arange(-40, 40, self.config.overlays['Contour']['Interval'])
            contours=ax_n.contour(lon, lat, contourdata, levels, linewidths=1, antialiased=True, alpha=1, colors='k',
                                  linestyles='solid', transform=self.Proj)
            if self.config.overlays['Contour']['Label']== True: 
                ax_n.clabel(contours, levels,  fmt='%r', inline=True, fontsize=10)

        # * * *
        ## PLOT BACKGROUND FEATURES
        if self.config.overlays['GEO']['OX'] == True:
            Data['P_Radar'].plot_topo(Master_Plt, ax_n)

        # * * *
        ## DEAL WITH COLORBARS
        # Attach colorbar to each subplot
        divider = make_axes_locatable(plt.gca())
        c_ax = divider.append_axes("bottom", "5%", pad="2%", axes_class=plt.Axes)
        sm = plt.cm.ScalarMappable(cmap=colormap, norm=matplotlib.colors.Normalize(vmin=vminb, vmax=vmaxb))
        sm._A = []
        cb = plt.colorbar(sm, cax=c_ax, orientation='horizontal', extend='both', drawedges=testing_testing, boundaries=testing_bound)
        cb.set_label(label=c_label, size=self.Radar_title_font['fontsize']) #, weight='bold')

        #set up the colorbar for the thermo colorline overlay
        c_ax2 = divider.append_axes("right", "5%", pad="2%", axes_class=plt.Axes)
        Thermo_cbar = plt.colorbar(TVARS.CS3, cax=c_ax2, label=TVARS.Tvar_lab, ticks=MaxNLocator(integer=True), use_gridspec=True)

        #  patches = []
        #  fancybox = mpatches.FancyBboxPatch((0,0),1.5,.5, boxstyle=mpatches.BoxStyle("Round", pad=0.02), color='red')
        #  patches.append(fancybox)
        #  collection = PatchCollection(patches)
        #  ax_n.add_collection(collection)
        #  box = AnchoredAuxTransformBox(ax_n.transData, loc='upper left')
        #  box = HPacker(children=[box1, box2], align="center", pad=0, sep=5)
        #  anchored_box = AnchoredOffsetbox(loc='lower left',child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),
                                 #  bbox_transform=ax.transAxes,borderpad=0.)
        #  ax_n.add_artist(box)
        #  ax_n.add_artist(anchored_box)

        # * * *
        ## SET UP LEGENDS
        if leg == True: #this means you are currently making the left subplot
            #add legend for platform markers
            l = ax_n.legend(loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(-0.07,.5),
                            markerscale=1.1, labelspacing=.7, handletextpad=.9, handlelength=.1)#, title="Platforms")
            l.set_title("Platforms", prop=self.leg_title_font)

        ###################
        if self.config.print_long == True: print('~~~~~~~~~~~Made it through radar_subplots~~~~~~~~~~~')


##########################################################################
###########################
## PLOTTING DEFINTIONS  ###
###########################
def plotting(config, Data, TVARS, start_comptime, tilt=None):
    print("*** Entering plotting() ")
    plotting_start_time = time.time()
    ''' Initial plotting defenition: sets up fig size, layout, font size etc and will call timeseries and radar subplots
    ----------
    INPUTS: Data: (dictionary)
                contains objects corresponding to all available datasets (radar, torus insitue platforms etc) and  which particular variable
                should be plotted on time series and colorlines
            print_long & e_test: (bool strings)
                True/False vars that control how much information is printed out to the terminal window. Helpful for debugging but can be overkill
            start_comptime:
                Time at which you first begain plotting this particular image (will help to report out how long it took to create image)
    '''
    if config.print_long == True: print('~~~~~~~~~~~Made it into plotting~~~~~~~~~~~~~~~~~~~~~')
    #initilize the plot object(will have info about plotlayout and such now)
    PLT = Master_Plt(Data, config)

    ## Make the Radar Subplots (if applicable)
    #### * * * * * * * * * * * * * * * * * ** * * * *
    if len(config.r_mom) != 0:
        if 'sim_vel' in config.r_mom:
            sd_test = singledop.SingleDoppler2D( L=30.0, radar=Data['P_Radar'].rfile, range_limits=[1, 6], filter_data=True,
                    filter_distance=5, thin_factor=[4, 12], sweep_number=Data['P_Radar'].swp, name_vr='velocity')

            display = singledop.AnalysisDisplay(sd_test)
            #  display = singledop.AnalysisDisplay(Data['P_Radar'].rfile.fields[sd_test])
            display.four_panel_plot(scale=400, legend=20, return_flag=False, thin=6,
                        levels=-30.0+2.0*np.arange(31), name_vr='vel_fix', name_dz='refl_fix', dz_cmap='pyart_HomeyerRainbow' ,
                                    xlim=[-21,21], ylim=[-21,21])
        else:
            mom_index=0
            for row in range(PLT.R_rows):
                for col in range(PLT.R_cols):
                    try:
                        mom=config.r_mom[mom_index]
                        make_rplt = True
                    except: make_rplt= False

                    if make_rplt == True:
                        print('Radar plot: '+ mom +', Outer GridSpec Pos: [0, :], SubGridSpec Pos: ['+str(row)+', '+str(col)+']')
                        ax_n= PLT.fig.add_subplot(PLT.r_gs[row, col], projection= PLT.R_Proj)
                        if row==0 and col==0: leg = True
                        else: leg = False
                        PLT.radar_subplots(mom, day, ax_n, PLT.fig, TVARS, Data, leg)
                        mom_index =mom_index+1
                    else: pass

    ## Make the Times Series Subplots (if applicable)
    #### * * * * * * * * * * * * * * * * * ** * * * *
    if len(config.Time_Series) != 0:
        for (subrow,), ts in np.ndenumerate(config.Time_Series):
            print('Time Series plot: '+ ts +', Outer GridSpec Pos: [1, :], SubGridSpec Pos: ['+str(subrow)+', :]')
            if subrow==0: ax_n= PLT.fig.add_subplot(PLT.ts_gs[subrow,:])
            else:  ax_n=PLT.fig.add_subplot(PLT.ts_gs[subrow,:], sharex=ax_n)
            PLT.time_series(ts, ax_n, PLT.fig, TVARS, Data)

    ## Finish plot
    #### * * * * *
    ### Get outdir and outname (aka file path and file name) for the correct image setup
    #if radar is included in the image
    if len(config.r_mom) != 0: 
        plt.suptitle(Data['P_Radar'].site_name+' '+str(tilt)+r'$^{\circ}$ PPI '+Data['P_Radar'].fancy_date_str, y=PLT.title_spacer, ha='center')
        #if both radar and timeseries are included in the image
        if len(config.Time_Series) != 0:
            file_string = '_'.join(config.Time_Series)
            if ('Thetae' in config.Time_Series) and ('Thetav' in config.Time_Series): plt_type = 'both_Thetas'
            elif (len(config.Time_Series) != 1) and ('Wind' in config.Time_Series): plt_type = config.overlays['Colorline']['Var']+'/Wind'
            else: plt_type = config.overlays['Colorline']['Var']
        #if only radar is included in the image (aka no platform info overlayed)
        else:
            file_string = 'r_only'
            if config.overlays['Contour']['Lines']== True:
                plt_type='Radar_only/ppi/contour'
            else:
                plt_type='Radar_only/ppi'


        #add the file name
        if Data['P_Radar'].name in pform_names('KA'):
            output_name = Data['P_Radar'].site_name+'_'+Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+'_'+file_string+'.png'
        else: #aka WSR88D or NOXP    
            output_name = Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+'_'+file_string+'_'+Data['P_Radar'].site_name+'.png'
        #path to file 
        rmom_string= '_'.join(config.r_mom)
        plt_dir=Data['P_Radar'].dir_name+'/'+plt_type+'/'+rmom_string+'/'+str(tilt )+'_deg/'+config.radar_controls['Centered_Pform']+'/'+str(
                    config.radar_controls['offsetkm'][config.Radar_Plot_Type])+'km/'

    #if only timeseries are included in the image
    elif len(config.Time_Series) != 0 and len(config.r_mom) == 0:
        plt.suptitle(day+' Time Series', y=title_spacer)
        
        #add the file name
        file_string = '_'.join(config.Time_Series)
        if config.ts_extent ==None: file_ender='_full.png'
        else: file_ender='.png'
        output_name= day+'_TS_'+file_string+file_ender
        #path to file 
        plt_dir='TSeries/'
    
    pts = [np.nan]
    if config.clicker == True:
        print('hello')

        plt.setp(plt.gca(), autoscale_on=False)
        print('hello2')
        plt.draw()
        print('hello3')
        #  plt.waitforbuttonpress()
        print('hello4')
        pts = plt.ginput(n=40,timeout=0) # look up ginput docs for some good guidance
        print('hello5')
        print(pts)
    
    ## * * *     
    #This is the directory path for the output file
    outdir = config.g_TORUS_directory+day+'/plots/'+plt_dir
    # Setup Function Cache for speedup, if path directory does not make the directory
    if not os.path.exists(outdir): Path(outdir).mkdir(parents=True, exist_ok=True)
    output_path_plus_name = outdir + output_name
    print(output_path_plus_name)

    #  PLT.fig.tight_layout()
    plt.savefig(output_path_plus_name, bbox_inches='tight', pad_inches=.3)
    timer(start_comptime, time.time())
    plt.close('all')   # close it
    PLT.fig.clf()          # no really... close it
    gc.collect()       # get tf outta here

    plotting_time = plotting_start_time - time.time()
    print('\a') ## Makes a ding noise
    if config.print_long == True: print('~~~~~~~~~~~made it through plotting~~~~~~~~~~~~~~~~~~')
    print('*** Exiting Plotting() after ', plotting_time, 'sec \n \n***************************************************************************************************')

##############################################################################################
def process_instruments(config, day):
    print('\nRead in Torus Platforms')
    #print(pform_names('ALL')) #list of ALL possible platforms
    subset_pnames = [] #array of platforms that actually have data for the time range we are plotting
    Data = {} #dictionary in which all the class objects will be stored (will contain platform data, locations, plotting variables etc)

    ## Read in the data for the TORUS Insitu platforms (if available)
    Data, subset_pnames = Add_to_DATA(config, day, 'TInsitu', Data, subset_pnames)

    print('\nRead in Stationary Platforms Arrays')
    ## Read in the data for the Stationary Array platforms (if available)
    Data, subset_pnames = Add_to_DATA(config, day, 'STN_I', Data, subset_pnames)
    return Data, subset_pnames


##############################################################################################
## Plot Radar ##
################
def plot_radar(config, day, Data, TVARS):
    print('\nYes Plot Radar \n')
    def find_radar_files(config, day, Data):
        r_files_path =[]
        if config.Radar_Plot_Type == 'KA_Plotting':
            ## Get radar files
            path = config.g_TORUS_directory + day+'/data/radar/TTUKa/netcdf/*/dealiased_*'
            r_files_path= sorted(glob.glob(path))

        elif config.Radar_Plot_Type == 'NOXP_Plotting':
            ## Get radar files
            path = config.g_TORUS_directory + day+'/data/radar/NOXP/'+day+'/*/sec/dealiased_ordered_*'
            #  path = config.g_TORUS_directory + day+'/data/radar/NOXP/'+day+'/*/sec/fixed*'
            r_files_path = sorted(glob.glob(path))

        elif config.Radar_Plot_Type == 'WSR_Plotting':
            r_files_path = det_nearest_WSR(Data[config.radar_controls['Centered_Pform']].df)
        return r_files_path
    
    ####
    def plot_radar_file(config, r_file, day, Data, TVARS, subset_pnames):
        print("\n*******************************\n plot_radar_file: scan file_name = {}\n".format(r_file))
        start_comptime = time.time()
        if config.Radar_Plot_Type == 'WSR_Plotting':
            is_WSR= True
            valid_time= True
        else:
            is_WSR= False 
            if config.Radar_Plot_Type == 'KA_Plotting':
                time_string= str(20)+str(os.path.split(r_file)[1][13:24])
                rtime = datetime.strptime(time_string, "%Y%m%d%H%M%S")
            elif config.Radar_Plot_Type == 'NOXP_Plotting':
                head_tail= os.path.split(r_file)
                rtime = datetime.strptime(head_tail[1][24:39], "%Y%m%d_%H%M%S")
            time_start, time_end = time_range(config)
            valid_time = time_in_range(time_start, time_end, rtime)

        # * * * 
        for tilt in config.radar_controls['p_tilt']:
            if valid_time == True:
                #open file using pyart
                radar = read_from_radar_file(r_file, is_WSR)
                swp_id = sweep_index(tilt, is_WSR, radar)
                if config.Radar_Plot_Type == 'KA_Plotting':
                    ## Only proceed if file is a ppi scan; we are currently not interested in plotting the RHI scans
                    if radar.scan_type != 'ppi': return
            # * * * 
            if (valid_time == False) or (swp_id == None):
                if valid_time == False:
                    print(r_file+' was not in the timerange being plotted')
                elif swp_id == None:
                    print(r_file+' did not have a sweep that matched our tilt angle of intrest ('+str(tilt)+')')
                #end evaluation for this file (do not make a plot)
            else:
                ## Read in radar data and add to Data dict
                ##### + + + + + + + + + + + + + + + + + + +
                #  Establish info for the main plotting radar (aka scantime etc) & and locations for other radars (if deployed)
                Data, subset_pnames = Add_to_DATA(config, day, 'RADAR', Data, subset_pnames, MR_file=radar, swp=swp_id)
                if config.print_long == True: print(str(Data)+'\n')

                ## Proceed to plot the radar
                ##### + + + + + + + + + + + +
                print("\nProducing Radar Plot:")
                plotting(config, Data, TVARS, start_comptime, tilt)
                print("done in plot_radar_file")


    #########
    if config.Radar_Plot_Type == 'WSR_Plotting':
        WSR_info = find_radar_files(config, day, Data)
        radar_files_paths= []
        for rad, times in WSR_info.items():
            print("Radar_Site: {}\n Start: {}\n End: {} \n***".format(rad, times[0], times[1]))
            radar_files = get_WSR_from_AWS(config, day, times[0], times[1], rad)
            #  radar_files_paths.append(radar_files)
            radar_files_paths= radar_files_paths+radar_files
    else:
        radar_files_paths = find_radar_files(config, day, Data)

    ## Proceed to plot the radar
    if config.print_radar_info== True: info_radar_files('compact', radar_files_paths )
    print('******\n Radar files to process:')
    pp.pprint(radar_files_paths)
    Parallel(n_jobs=config.nCPU, verbose=10)(delayed(plot_radar_file)(config, r_file, day, Data, TVARS, subset_pnames) for r_file in radar_files_paths)


## Plot Time Series##
#####################
def plot_time_series(config, day, Data, TVARS):
    print("\nPlot Timeseries\n")
    start_comptime = time.time()
    plotting(config, day, Data, TVARS, start_comptime)

def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

##################################################
# Define a triangle by clicking three points


plt.clf()
plt.setp(plt.gca(), autoscale_on=False)

tellme('You will define a triangle, click to begin')

plt.waitforbuttonpress()

while True:
    pts = []
    while len(pts) < 3:
        tellme('Select 3 corners with mouse')
        pts = np.asarray(plt.ginput(3, timeout=-1))
        if len(pts) < 3:
            tellme('Too few points, starting over')
            time.sleep(1)  # Wait a second

    ph = plt.fill(pts[:, 0], pts[:, 1], 'r', lw=2)

    tellme('Happy? Key click for yes, mouse click for no')

    if plt.waitforbuttonpress():
        break

    # Get rid of fill
    for p in ph:
        p.remove()


##################################################
# Now contour according to distance from triangle
# corners - just an example

# Define a nice function of distance from individual pts
def f(x, y, pts):
    z = np.zeros_like(x)
    for p in pts:
        z = z + 1/(np.sqrt((x - p[0])**2 + (y - p[1])**2))
    return 1/z


X, Y = np.meshgrid(np.linspace(-1, 1, 51), np.linspace(-1, 1, 51))
Z = f(X, Y, pts)

CS = plt.contour(X, Y, Z, 20)

tellme('Use mouse to select contour label locations, middle button to finish')
CL = plt.clabel(CS, manual=True)

##################################################
# Now do a zoom

tellme('Now do a nested zoom, click to begin')
plt.waitforbuttonpress()

while True:
    tellme('Select two corners of zoom, middle mouse button to finish')
    pts = plt.ginput(2, timeout=-1)
    if len(pts) < 2:
        break
    (x0, y0), (x1, y1) = pts
    xmin, xmax = sorted([x0, x1])
    ymin, ymax = sorted([y0, y1])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

tellme('All Done!')
plt.show()

#############################################################################################
for day in plot_config.day_list:
    #####################################################
    # Read in Data that does not change for each image ##
    #####################################################
    Data, subset_pnames = process_instruments(plot_config, day)
    ## Establish info for the plotting variable (aka max min etc)
    TVARS=Thermo_Plt_Vars(plot_config, day, Data)

    ###################################
    # Create Plot for each radarfile ##
    ###################################
    if len(plot_config.r_mom) != 0:
        plot_radar(plot_config, day, Data, TVARS)

    ################################
    # Create Timeseries only plot ##
    ################################
    elif len(plot_config.Time_Series) != 0:
        # Create Timeseries only plot ##
        #Only plot timeseries (this code isn't fully fleshed out but in theroy this code is built in such a way to allow for this)
        plot_time_series(plot_config, day, Data, TVARS)

###########################
plt.rcdefaults()
timer(totalcompT_start, time.time(), total_runtime=True)
print("ALL FINISHED")

tracemalloc.stop()

var_list = ['T2', 'T3']
Num_of_members=2
Data={}
for var in var_list:
    ES_top=0
    for mem_num in arange(1,Num_of_members):
        file_name='filename...'+str(mem_nm)

        ds =np.Dataset(file_name)
        fcst=ds[var][:]
        if ES_top ==0:
            ES_top

        ES_top=ES_top+fcst

    ES_mean = ES_top/Num_of_members
    Data.update({var:ES_mean})

print(Data['T2'])

def example(species, real_name, asker='Ellie'):
    if species=='dog':
        name= real_name
    elif species == 'cat':
       if asker=='Ellie': 
           name= 'Humphrey'
       else:
           name= real_name
    else: 
        name= 'Try again'
    return name 
real_n='Chester'
A= 'Jillian'
what_called= example('dog', real_n)
what_called_2= example('dog', real_n, A)
print(what_called)

