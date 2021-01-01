#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import needed modules
######################
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
#  import matplotlib.font_manager as mfm
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
#  import modin.pandas as pd
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

import tracemalloc
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
if plot_config.country_roads == True: import osmnx as ox

## Read in defns that I have stored in another file (for ease of use/consistancy accross multiple scripts)
from read_pforms import pform_names, Add_to_DATA, Platform, Torus_Insitu, Radar, Stationary_Insitu
from read_pforms import error_printing, timer, time_in_range, time_range

totalcompT_start = time.time()
################################################################################################
##########
# Classes
##########
class Thermo_Plt_Vars:
    def __init__ (self, config, Data):
        self.Te_lab, self.Tv_lab, self.Tf_lab = "Equi. Pot. Temp [K]", "Vir. Pot. Temp [K]", "Temp (fast) [K]"
        self.Tv_GMin, self.Tv_GMax = self.find_global_max_min('Thetav', Data)
        self.Te_GMin, self.Te_GMax = self.find_global_max_min('Thetae', Data)
        self.Tf_GMin, self.Tf_GMax = self.find_global_max_min('tfast', Data)

        #settings for the thermo var being plotted on radar subplots
        if len(config.r_mom) != 0:
            if config.R_Tvar == "Thetae":
                lab, GMin, GMax = self.Te_lab, self.Te_GMin, self.Te_GMax
            elif config.R_Tvar == "Thetav":
                lab, GMin, GMax = self.Tv_lab, self.Tv_GMin, self.Tv_GMax
            elif config.R_Tvar == "tfast":
                lab, GMin, GMax = self.Tf_lab, self.Tf_GMin, self.Tf_GMax
            self.Tvar_lab, self.Tvar_GMin, self.Tvar_GMax = lab, GMin, GMax

            ## Make a dummy plot to allow for ploting of colorbar
            #  cmap=plt.get_cmap('rainbow')
            self.norm =plt.Normalize(GMin, GMax)

            Z = [[0,0],[0,0]]
            levels = np.arange(GMin, GMax+1, 1)

            self.CS3 = plt.contourf(Z, levels, cmap=cmocean.cm.thermal)
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
        #  print(global_max)
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
        plt.rc('font', size= 15)         # controls default text sizes
        plt.rc('axes', labelsize= 21) # fontsize of the axes title, and x and y labels
        plt.rc('legend', fontsize= 23, borderpad=.5, facecolor='white', edgecolor= 'black', shadow=True, fancybox=True, framealpha=1)# legend fontsize
        plt.rc('figure', titlesize= 50, facecolor='white')  # fontsize of the figure title
        #  self.leg_title_font=FontProperties(size=25, weight='bold')
        self.leg_title_font={'size':25, 'weight':'bold'}
        self.Radar_title_font= {'fontsize':40}
        #  plt.rc('savefig', bbox= 'tight', pad_inches=.3)
        #  plt.rc('lines', markeredgewidth= 'grey')
        #  plt.rc('path.effects', )

        ## Establish the Plot size and layout
        #### * * * * * * * * * * * * * * * * *
        # This is the layout for radar and time series subplots
        if len(config.Time_Series) != 0 and len(config.r_mom) != 0:
            self.fig= plt.figure(figsize=(32,20))
            if len(config.Time_Series) == 1:
                self.outer_gs= GridSpec(nrows=2, ncols=1, height_ratios=[3,2], hspace=.1)
            if len(config.Time_Series) == 2:
                self.outer_gs= GridSpec(nrows=2, ncols=1, height_ratios=[4,3], hspace=.1)
                if 'Wind' in config.Time_Series: hratio=[2,3]
                else: hratio=[1,1]
            self.ts_gs = GridSpecFromSubplotSpec(len(config.Time_Series), 1, subplot_spec=self.outer_gs[1, :], height_ratios=hratio, hspace=0)
            self.r_gs = GridSpecFromSubplotSpec(1, len(config.r_mom), subplot_spec=self.outer_gs[0, :])

        # This is the layout for time series only
        if len(config.Time_Series) != 0 and len(config.r_mom) == 0:
            print('this is the layout for time series only')
            #  self.fig= plt.figure(figsize=(32,5))
            self.fig= plt.figure(figsize=(32,10))
            if len(config.Time_Series) == 1:  hratio=[1]
            if len(config.Time_Series) == 2:
                if 'Wind' in config.Time_Series: hratio=[2,3]
                else: hratio=[1,1]
            self.ts_gs = GridSpec(len(config.Time_Series), 1, height_ratios=hratio, hspace=0)

        # This is the layout for radar only
        if len(config.r_mom) != 0 and len(config.Time_Series) == 0:
            print('this is the layout for radar plots only')

        ## Establish some vars that will be helpful to the plot later on
        #### * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        #if we are plotting radar
        if len(config.r_mom) != 0:
            # The extent of the area to be plotted
            #redefine the classvariable for Domain and display
            self.Domain = Platform.getLocation(Data[config.Centered_Pform], scan_time=Data['P_Radar'].Scan_time, offsetkm= config.offsetkm)
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
            if self.config.ts_extent != None:
                ax.set_xlim(scan_time - timedelta(minutes=self.config.ts_extent), scan_time + timedelta(minutes=self.config.ts_extent))

            # Yaxis
            if ts in ['Thetav','Thetae']:
                if self.config.man_ylim == True:
                    ax.set_ylim(335,350)
                else:     
                    ax.set_ylim(TVARS.Tvar_GMin, TVARS.Tvar_GMax)

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
            ax.xaxis.set_minor_locator(AutoMinorLocator(6))

            # Yaxis
            if ts in ['Thetav','Thetae']:
                ax.yaxis.set_major_locator(MultipleLocator(5)) # set up major tick marks (this is set up to go by 5's will want to change for diff vars)
                ax.yaxis.set_minor_locator(AutoMinorLocator(5)) # set up minor ticks (this have it increment by 1's will want to change for diff vars)
            elif ts == 'Wind' and twinax == False:  ax.yaxis.set_major_locator(LinearLocator(numticks=5))
            elif ts == 'Wind' and twinax == True:  ax.yaxis.set_major_locator(FixedLocator(np.arange(0, 450, 90)))

        #if you are calling this for a radar subplot locate the ticks
        if radar != None and interval != None:
            def calc_ticks(a,b,interval):
                ''' Given: range from a to b where a is negative and b is positive
                    Return: an array of tick locations including 0 and spaced interval amounts above and below 0
                '''
                left, right = np.arange(0, a-1, -interval), np.arange(0, b+1, interval)
                ticks = np.concatenate([left,right])
                ticks = np.sort(ticks)
                ticks = np.unique(ticks)
                return ticks
            #Det the boundaries of the domain
            x0, x1 = ax.get_xbound()
            y0, y1 = ax.get_ybound()

            #Det the tick locations
            yticks, xticks = calc_ticks(y0, y1, interval), calc_ticks(x0, x1, interval)
            #set the calculated ticks on the graph
            ax.set_yticks(yticks)
            ax.set_xticks(xticks)
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
            ax.tick_params(which='both', axis='both', direction='in', top=True, bottom=False, labelbottom=False,labeltop=True)
            ax.tick_params(which='both', axis='y', labelrotation=90)
        # + + + + + + + + + + + + ++ + +

        ##Grid Display Characteristics
        # # # # # # # # # # # # # # # #
        if ts != None:
            ax.tick_params(which='major', axis='both', grid_linewidth=2.5)
            ax.tick_params(which='major', axis='y', grid_color='grey', grid_linewidth=2, grid_linestyle='--', grid_alpha=.8)

        if radar != None:
            ax.tick_params(which='both', axis='both', grid_linewidth=2, grid_color='grey', grid_alpha=.8)
        ax.tick_params(which='minor', axis='both', grid_linestyle=':', grid_alpha=.8)
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
                return '{} km'.format(x)
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
        if ts in ['Thetav', 'Thetae']:
            for p in Data.values():
                if isinstance(p, Torus_Insitu):
                    if self.config.print_long == True: print('Plotting '+str(p.name)+' on time series')
                    ## If the platform matches a type listed in TS_masking use masked data for the time series; else use the unmasked data
                    if p.type in self.config.TS_masking[:]:
                        if ts == "Thetae":   plotting_data= p.Thetae_ts_mask_df
                        elif ts == "Thetav":   plotting_data= p.Thetav_ts_mask_df
                    else: plotting_data= p.df[ts].values

                    ## Plot
                    ax_n.plot(p.df['datetime'], plotting_data, linewidth=3, color=p.l_color, label=p.leg_str) #assigning label= is what allows the legend to work

                    TSleg_entry = Line2D([], [], label=p.leg_str, linewidth=12, color=p.l_color)
                    TSleg_elements.append(TSleg_entry)
            ## Set up XY axes tick locations and Labels
            self.tick_grid_settings(ax=ax_n, ts=ts, scan_time=Data['P_Radar'].Scan_time)
            
            if ts == "Thetae":   ylab = TVARS.Te_lab
            elif ts == "Thetav":   ylab = TVARS.Tv_lab
            ax_n.set_ylabel(ylab)


        # * * *
        if ts == 'Wind':
            p = Data[self.config.Wind_Pform]
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
        #  leg = ax_n.legend(handles= TSleg_elements, loc='center left')
        leg = ax_n.legend(handles= TSleg_elements, scatterpoints=3, loc='center left')
        if ts == 'Wind':
            leg.set_title(self.config.Wind_Pform, prop=self.leg_title_font)
            leg.remove()
            ax_2.add_artist(leg)

        #if plotting more than one time series then only include the x axis label and ticks to the bottom timeseries
        #only label ticks for the outer axis of the timeseries (aka remove overlapping axis)
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
            ax_n.axvspan(Data['P_Radar'].Scan_time - timedelta(minutes=self.config.cline_extent), 
                         Data['P_Radar'].Scan_time + timedelta(minutes=self.config.cline_extent), facecolor='0.5', alpha=0.4)
        ##Label Formatter
        # # # # # # # # #
        ax_n.set_xlabel('Time (UTC)')
        ax_n.yaxis.set_label_coords(-.03, .5)
        '''
        print('TS ', ts)
        bbox = ax_n.get_tightbbox(fig.canvas.renderer, call_axes_locator=True)
        x0, y0, width, height = bbox.transformed(fig.transFigure.inverted()).bounds
        # slightly increase the very tight bounds:
        #  xpad, ypad = width, height 
        #  ax_n.patch((x0, y0), width, height, edgecolor='red', linewidth=3, fill=False)
        #  box = ax_n.get_position()
        #  box2=ax_n.get_tightbbox(fig.canvas.get_renderer(), call_axes_locator = True)
        '''
        if self.config.print_long == True: print('~~~~~~~~~~~Made it through time_series~~~~~~~~~~~~~~')


    # * * * * * *  *
    def radar_subplots(self, mom, ax_n, fig, TVARS, Data, leg):
        ''' Plots each of the radar subplots including the marker overlays corresponding to the additional platforms
        ----
        INPUTS:
        mom: str, indicates which radar moment you would like to plot (ie Reflectivity, Velocity, Spectrum Width etc)
        Data: dict as described in the plotting defn
        leg: bool str whether or not you want a legend associated with this particular subplot
        print_long & e_test: bool str as described in the ppi defn
        '''
        if self.config.print_long == True: print('~~~~~~~~~~~made it into radar_subplots~~~~~~~~~~~~~~')

        ## SET UP VARS FOR EACH RADAR MOMENTS
        if mom == 'refl':
            p_title, c_label, c_scale = 'Reflectivity', 'Radar Reflectivity [dbz]', 'pyart_HomeyerRainbow'
            if Data['P_Radar'].name in pform_names('KA'):
                field, vminb, vmaxb, sweep = 'refl_fix', -30., 30., Data['P_Radar'].swp
            if Data['P_Radar'].name == 'NOXP':
                field, vminb, vmaxb, sweep = 'DBZ', -10., 60., 0# Data['P_Radar'].swp
            if Data['P_Radar'].name == 'WSR88D':
                field, vminb, vmaxb, sweep =  'reflectivity', -10., 75., Data['P_Radar'].swp[0]
        elif mom == 'vel':
            p_title, c_label, c_scale = 'Radial Velocity', 'Velocity [m/s]', 'pyart_balance'
            if Data['P_Radar'].name in pform_names('KA'):
                field, vminb, vmaxb, sweep = 'vel_fix', -30., 30., Data['P_Radar'].swp
            if Data['P_Radar'].name == 'NOXP':
                field, vminb, vmaxb, sweep = 'VEL', -40., 40., 0  #Data['P_Radar'].swp
            if Data['P_Radar'].name == 'WSR88D':
                field, vminb, vmaxb, sweep = 'velocity', -40., 40., Data['P_Radar'].swp[1]
        else: print("Hey what just happened!\n Check the Radar moments for spelling")

        ## Plot the radar
        ax_n.set_title(p_title, y=-.067, fontdict=self.Radar_title_font)
        self.display.plot_ppi_map(field, sweep, ax=ax_n, cmap=c_scale, vmin=vminb, vmax=vmaxb, width=self.config.offsetkm*2000, 
                                  height=self.config.offsetkm*2000, title_flag=False, colorbar_flag=False, embelish=False)

        # Has to be here or it doesn't work
        ax_n.set_extent(self.Domain)
        self.tick_grid_settings(ax=ax_n, radar=True,scan_time=None, interval=10*1000)
        ax_n.grid(True)
        # Do it again here to restore extent
        ax_n.set_extent(self.Domain)

        # * * *
        ## PLOT PLATFORMS AS OVERLAYS(ie marker,colorline etc) ON RADAR
        #  iterate over each object contained in dict Data (returns the actual objects not their keys)
        def plot_markers(p, lon, lat, lab_offset=(0,0), lab_c='xkcd:pale blue', zipper_lab=None):
            ax_n.plot(lon, lat, transform=self.Proj, marker=p.m_style, mfc=p.m_color, linestyle=None, 
                      mew=2, ms=p.m_size, mec='k',label=p.leg_str, zorder=10)
            
            ## Optional textlabel on plot
            if p.marker_label == True:
                if zipper_lab == None:
                    ax_n.text(p.lon+lab_offset[0], p.lat-lab_offset[1], p.name, transform=self.Proj,
                              path_effects=[PE.withStroke(linewidth=4, foreground=lab_c)])
                else:
                    for x, y, lab in zip(p.lon, p.lat, zipper_lab):
                        trans = self.R_Proj.transform_point(x+.009, y-.002, self.Proj)
                        ax_n.text(trans[0], trans[1], lab, fontsize= 20)

        def plot_TORUSpform(Data, p, border_c='xkcd:light grey', labelbias=(0,0)):
            #grab the subset of data of +- interval around radar scan
            p_sub, p_deploy = p.grab_pform_subset(p, Data, time_offset=self.config.cline_extent)

            #if there is data for the platform that falls within the time and location of interest
            if p_deploy == True:
                ##READ in data from p_sub 
                x, y, C = p_sub['lon'].values, p_sub['lat'].values, p_sub[self.config.R_Tvar].values

                ##PLOT the line that extends +/- min from the platform location; 
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

                #  determine the row of the pform data at scantime (for plotting marker)
                C_Point=p_sub.loc[p_sub['datetime'].sub(Data['P_Radar'].Scan_time).abs().idxmin()]
                return C_Point.lon, C_Point.lat 

            elif p_deploy == False:
                if self.config.print_long == True: print('The platform was not deployed at this time')
                return False, False

        def rhi_spokes_rings(pform, ax= ax_n):
            ''' Plot the RHI spoke and ring for a radar
            '''
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
        for p in Data.values():
            if self.config.print_long == True: print(p.name)
            ######
            #Plot Inistu Torus Platforms (if desired and available)
            if isinstance(p, Torus_Insitu):
                X, Y = plot_TORUSpform(Data, p, border_c='k')
                if X != False:
                    ##PLOT the platform marker 
                    plot_markers(p, X, Y)

            ######
            #Plot Stationary Inistu Platforms ( aka mesonets and ASOS) if desired and available
            if isinstance(p, Stationary_Insitu):
                # determine if any of the sites fall within the plotting domain
                sites_subdf, valid_sites = p.grab_pform_subset(p, Data, bounding= self.Domain)
                # if there are sites within the domain plot the markers and include in the legend
                if valid_sites == True:
                    plot_markers(p, sites_subdf.lon, sites_subdf.lat, lab_offset=(.009,.002), zipper_lab = sites_subdf.Stn_ID)

            #####
            #Plot Radar Platforms (if desired and available)
            if isinstance(p, Radar):
                if p.type != 'MAINR':
                    #det if (any of) the radar is located within the area included in the plot domain at the time of the plot
                    sites_subdf, valid_sites = p.grab_pform_subset(p, Data, bounding= self.Domain)

                    #if these conditions are met then plot the radar marker(s)
                    if valid_sites == True:
                        if p.type == 'KA' or p.type == 'NOXP':
                            ## Plot the marker
                            plot_markers(p, p.lon, p.lat, lab_offset=(.009,.002))

                            ## Plot RHI spokes
                            if p.type == 'KA':
                                if self.config.rhi_ring == True: rhi_spokes_rings(p)

                        if p.type == 'WSR':
                            ## Plot the marker
                            plot_markers(p, sites_subdf.lon, sites_subdf.lat, lab_offset=(.03,.01), zipper_lab=sites_subdf)

        
        ## PLOT BACKGROUND FEATURES
        Data['P_Radar'].plot_topo(Master_Plt, ax_n)
        ## DEAL WITH COLORBARS
        # Attach colorbar to each subplot
        divider = make_axes_locatable(plt.gca())
        c_ax = divider.append_axes("right", "5%", pad="2%", axes_class=plt.Axes)
        sm = plt.cm.ScalarMappable(cmap=c_scale, norm=matplotlib.colors.Normalize(vmin=vminb, vmax=vmaxb))
        sm._A = []
        cb = plt.colorbar(sm, cax=c_ax, label=c_label)

        ## SET UP LEGENDS
        if leg == True: #this means you are currently making the left subplot
            #add legend for platform markers
            #  l = ax_n.legend(handles=legend_elements, loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(-0.1,.5), handlelength=.1)#, title="Platforms")
            l = ax_n.legend(loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(-0.07,.5), markerscale=1.2, labelspacing=.6, handlelength=.1)#, title="Platforms")
            #  h, le = ax_n.get_legend_handles_labels()
            #  l = ax_n.legend(loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(0,.5), handlelength=.1)#, title="Platforms")
            l.set_title("Platforms", prop=self.leg_title_font)
            #  print(vars(l))
            #  lbbox=l.get_window_extent()
            #  print(lbbox)
        if leg == False:  #this means you are currently making the right subplot
            ## Plot platform colorbar
            #set up colorbar axis that will be as tall and 5% as wide as the 'parent' radar subplot
            cbar_ax = inset_axes(ax_n, width= '5%', height= '100%', loc='center left', bbox_transform=ax_n.transAxes, bbox_to_anchor=(-.226, 0,1,1))
            cbar = plt.colorbar(TVARS.CS3, cax=cbar_ax, orientation='vertical', label=TVARS.Tvar_lab, ticks=MaxNLocator(integer=True))#,ticks=np.arange(Data['p_var'].global_min, Data['p_var'].global_max+1,2))
        '''
        print('RMOM ', mom)
        box = ax_n.get_position()
        box2=ax_n.get_tightbbox(fig.canvas.get_renderer(), call_axes_locator = True)
        ax_n.set_frame_on(True)
        '''
        # Shrink current axis by 20%
        #  box = ax.get_position()
        #  ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        # Shrink current axis's height by 10% on the bottom
        #  box = ax.get_position()
        #  ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
          #  plt.subplots_adjust(right=0.7)
        ###################
        if self.config.print_long == True: print('~~~~~~~~~~~Made it through radar_subplots~~~~~~~~~~~')


##########################################################################
###########################
## PLOTTING DEFINTIONS  ###
###########################
def plotting(config, Data, TVARS, start_comptime):
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
        for (subcol,), mom in np.ndenumerate(config.r_mom):
            print('Radar plot: '+ mom +', Outer GridSpec Pos: [0, :], SubGridSpec Pos: [:, '+str(subcol)+']')
            ax_n= PLT.fig.add_subplot(PLT.r_gs[:, subcol], projection= PLT.R_Proj)
            print(ax_n)
            if subcol==0: leg = True
            else: leg = False
            PLT.radar_subplots(mom, ax_n, PLT.fig, TVARS, Data, leg)

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
    #Add a title in the plot
    title_spacer=.93

    if len(config.Time_Series) != 0:
        #PLT TITLE STUFF
        file_string = '_'.join(config.Time_Series)
        if config.ts_extent !=None:
            output_name= config.day+'_TS_'+file_string+'.png'
        if config.ts_extent ==None:
            output_name= config.day+'_TS_'+file_string+'_full.png'

        # OUTDIR name stuff
        if ('Thetae' in config.Time_Series) and ('Thetav' in config.Time_Series): Thermo_type='both'
        else: Thermo_type= config.R_Tvar
        
        if len(config.r_mom) != 0:
            #PLT TITLE STUFF
            plt.suptitle(config.day+' Time Series', y=title_spacer)
            plt_dir=Data['P_Radar'].dir_name+'/'+Thermo_type+'/'+str(config.p_tilt)+'_deg/'
        else:
            plt_dir='TSeries/'+Thermo_type+'/'
        
        #This is the directory path for the output file
        outdir = config.g_plots_directory+config.day+'/plots/'+plt_dir+config.Wind_Pform+'/'
        # Setup Function Cache for speedup, if path directory does not make the directory
        if not os.path.exists(outdir): Path(outdir).mkdir(parents=True)
    
    if len(config.r_mom) != 0:
        if Data['P_Radar'].name in pform_names('KA') or Data['P_Radar'].name =='WSR88D':
            plt.suptitle(Data['P_Radar'].site_name+' '+str(config.p_tilt)+r'$^{\circ}$ PPI '+Data['P_Radar'].fancy_date_str, y=title_spacer)
        else:
            plt.suptitle(Data['P_Radar'].name+' '+str(config.p_tilt)+r'$^{\circ}$ PPI '+Data['P_Radar'].fancy_date_str, y=title_spacer)
        
        #add the file name
        if Data['P_Radar'].name in pform_names('KA'):
            output_name= Data['P_Radar'].site_name+'_'+Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+'_'+file_string+'.png'
        if Data['P_Radar'].name == 'NOXP':
            output_name = Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+'_'+file_string+'_NOXP.png'
        if Data['P_Radar'].name == 'WSR88D':
            output_name = Data['P_Radar'].Scan_time.strftime('%m%d_%H%M')+'_'+file_string+'_'+Data['P_Radar'].site_name+'.png'
    
    output_path_plus_name = outdir+output_name
    print(output_path_plus_name)

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
#######################
## RADAR DEFINITIONS ##
#######################
def read_from_radar_file(radar_file, Rtype):
    if Rtype== 'WSR':  radar = pyart.io.read_nexrad_archive(radar_file)
    elif Rtype== 'KA': radar = pyart.io.read(radar_file)
    elif Rtype =='NOXP': radar = pyart.io.read(radar_file)
    return radar
function_cache_memory = Memory(plot_config.g_cache_directory,verbose=1)
cached_read_from_radar_file = function_cache_memory.cache(read_from_radar_file)

# * * *
def plot_radar_file(config, r_file, Data, TVARS, subset_pnames):
    def sweep_index(radar= None):
        swp_id = None
        if config.Radar_Plot_Type == 'WSR_Plotting':
            #Hard code the swp numbers that will be associated with a given tilt angle
            if config.p_tilt == .5: swp_id=[0 , 1]
            elif config.p_tilt == 1: swp_id=[2 , 3]
            elif config.p_tilt == 1.5: swp_id=[4 , 5]
            else: print('The tilt angle {} is not hard coded yet for WSR'.format(config.p_tilt))

        else: 
            for i in range(radar.nsweeps):
                tilt_ang = radar.get_elevation(i) ## Det the actual tilt angle of a given sweep (returns an array)
                #  print("tilt: ", tilt_ang[0], ' p_tilt: ', config.p_tilt)
                #  print(np.around(tilt_ang[0], decimals=1))
                ## Check to see if the radarfile matches the elevation tilt we are interested in
                if np.around(tilt_ang[0], decimals=1) == config.p_tilt: swp_id = i 
        return swp_id
    #####
    print("\n*******************************\n plot_radar_file: scan file_name = {}\n".format(r_file))
    start_comptime = time.time()
    
    #####
    if config.Radar_Plot_Type == 'WSR_Plotting': 
        rad_type = 'WSR'
        valid_time= True
    
    elif config.Radar_Plot_Type == 'KA_Plotting':
        rad_type = 'KA'
        time_string= str(20)+str(os.path.split(r_file)[1][13:24])
        rtime = datetime.strptime(time_string, "%Y%m%d%H%M%S")
        time_start, time_end = time_range(config)
        valid_time = time_in_range(time_start, time_end, rtime)
        if valid_time == True:
            ## If file contains a ppi scan proceed; we are currently not interested in plotting the RHI scans
            if radar.scan_type == 'ppi':
                valid_time=False

    elif config.Radar_Plot_Type == 'NOXP_Plotting':
        rad_type = 'NOXP'
        head_tail= os.path.split(r_file)
        rtime = datetime.strptime(head_tail[1][6:21], "%Y%m%d_%H%M%S")
        time_start, time_end = time_range(config)
        valid_time = time_in_range(time_start, time_end, rtime)
        
    #####
    if valid_time == True:
        #open file using pyart
        radar = cached_read_from_radar_file(r_file, rad_type)
        swp_id = sweep_index(radar)      
    #####
    if (valid_time == False) or (swp_id == None):
        if valid_time == False:        
            print(r_file+' was not in the timerange being plotted')
        if swp_id == None:
            print(r_file+' did not have a sweep that matched our tilt angle of intrest')
        #end evaluation for this file (do not make a plot)
        return

    else:
        ## Read in radar data and add to Data dict
        ##### + + + + + + + + + + + + + + + + + + +
        #  Establish info for the main plotting radar (aka scantime etc) & and locations for other radars (if deployed)
        Data, subset_pnames = Add_to_DATA(config, 'RADAR', Data, subset_pnames, MR_file=radar, swp=swp_id)
        if config.print_long == True: print(str(Data)+'\n')
        if config.print_radar_info== True: print(radar.info(level='compact'))

        ## Proceed to plot the radar
        ##### + + + + + + + + + + + +
        print("\nProducing Radar Plot:")
        plotting(config, Data, TVARS, start_comptime)
        print("done in plot_radar_file")

##############################################################################################
def process_instruments(config):
    print('\nRead in Torus Platforms')
    #print(pform_names('ALL')) #list of ALL possible platforms
    subset_pnames = [] #array of platforms that actually have data for the time range we are plotting
    Data = {} #dictionary in which all the class objects will be stored (will contain platform data, locations, plotting variables etc)

    ## Read in the data for the TORUS Insitu platforms (if available)
    Data, subset_pnames = Add_to_DATA(config, 'TInsitu', Data, subset_pnames)

    print('\nRead in Stationary Platforms Arrays')
    ## Read in the data for the Stationary Array platforms (if available)
    Data, subset_pnames = Add_to_DATA(config, 'STN_I', Data, subset_pnames)
    return Data, subset_pnames

def find_radar_files(config, Data):
    r_files_path =[]
    # * * *
    if config.Radar_Plot_Type == 'KA_Plotting':
        ## Get radar files
        path = config.g_mesonet_directory + config.day+'/radar/TTUKa/netcdf/*/dealiased_*'
        r_files_path= sorted(glob.glob(path))
        if len(radar_files)==0:
            path = '/Users/severe2/Research/TORUS_Data/' +config.day+'/radar/TTUKa/netcdf/*/dealiased_*'
            r_files_path = sorted(glob.glob(path))
            if len(r_files_path)==0: print(path)

    # * * *
    elif config.Radar_Plot_Type == 'NOXP_Plotting':
        ## Get radar files
        path = config.g_mesonet_directory + config.day+'/radar/NOXP/'+config.day+'/*/sec/*'
        r_files_path = sorted(glob.glob(path))

    # * * *
    elif config.Radar_Plot_Type == 'WSR_Plotting':
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
        r_files_path = det_nearest_WSR(Data[config.Centered_Pform].df)
    return r_files_path

##############################################################################################
def plot_radar_files(config, Data, TVARS):
    def get_WSR_from_AWS(config, start, end, radar_id):
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
        path =  config.g_download_directory+ config.day +'/radar/Nexrad/Nexrad_files/'
        # If you dont have the path already make it and download the files
        if not os.path.exists(path): Path(path).mkdir(parents=True)
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

    print('\nYes Plot Radar \n')
    if config.Radar_Plot_Type == 'WSR_Plotting':
        WSR_info = find_radar_files(config, Data)
        for rad, times in WSR_info.items():
            print("Radar_Site: {}\n Start: {}\n End: {} \n***".format(rad, times[0], times[1]))
            radar_files = get_WSR_from_AWS(config, times[0], times[1], rad)
            ## Proceed to plot the radar
            print('********\n Radar files to process:')
            pp.pprint(radar_files)
            Parallel(n_jobs=config.nCPU, verbose=10)(delayed(plot_radar_file)(config, r_file, Data, TVARS, subset_pnames) for r_file in radar_files)

    else:
        radar_files_paths = find_radar_files(config, Data)
        ## Proceed to plot the radar
        print('******\n Radar files to process:')
        pp.pprint(radar_files_paths)
        Parallel(n_jobs=config.nCPU, verbose=10)(delayed(plot_radar_file)(config, r_file, Data, TVARS, subset_pnames) for r_file in radar_files_paths)

def plot_time_series(config, Data, TVARS):
    print("\nPlot Timeseries\n")
    #  + str(config.g_mesonet_directory+config.day+'/mesonets/NSSL/*.nc'))
    start_comptime = time.time()
    plotting(config, Data, TVARS, start_comptime)


#####################################################
# Read in Data that does not change for each image ##
#####################################################
Data, subset_pnames = process_instruments(plot_config)
## Establish info for the plotting variable (aka max min etc)
TVARS=Thermo_Plt_Vars(plot_config, Data)

###################################
# Create Plot for each radarfile ##
###################################
if len(plot_config.r_mom) != 0:
    plot_radar_files(plot_config, Data, TVARS)

################################
# Create Timeseries only plot ##
################################
elif len(plot_config.Time_Series) != 0:
    # Create Timeseries only plot ##
    #Only plot timeseries (this code isn't fully fleshed out but in theroy this code is built in such a way to allow for this)
    plot_time_series(plot_config, Data, TVARS)
###########################
plt.rcdefaults()
timer(totalcompT_start, time.time(), total_runtime=True)
print("ALL FINISHED")

tracemalloc.stop()
