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

register_matplotlib_converters()
logging.basicConfig(filename='this_is.log')
log = logging.getLogger(__name__)
# To run with a) warnings and b) stack trace on abort
# python3 -Walways  -q -X faulthandler plot_nexrad_insitu.py

## Imports form other files
############################
import config #this is the file with the plotting controls to access any of the vars in that file use config.var
if config.country_roads == True: import osmnx as ox

#rename a few commonly used vars so that the config.var does not have to be used repeatedly
print_long, e_test = config.print_long, config.e_test

## Read in defns that I have stored in another file (for ease of use/consistancy accross multiple scripts)
from read_pforms import Add_to_DATA, time_in_range, pform_names, Platform, Torus_Insitu, Radar, Stationary_Insitu, Pvar, error_printing, timer

totalcompT_start = time.time()
################################################################################################
def scale_bar(ax, length=None, location=(.96, -0.05), linewidth=3):
    """ ax is the axes to draw the scalebar on.
    length is the length of the scalebar in km.
    location is center of the scalebar in axis coordinates.
    (ie. 0.5 is the middle of the plot)
    linewidth is the thickness of the scalebar. """
    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())

    #Make tmc horizontally centred on the middle of the map,
    #vertically at scale bar location
    sbllx = (llx1 + llx0) / 2
    sblly = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(sbllx, sblly)

    #Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)

    #Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length:
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length)

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    #Plot the scalebar
    ax.plot(bar_xs, [sby, sby], transform=tmc, color='black', linewidth=linewidth, zorder=7)
    #Plot the scalebar label
    ax.text(sbx, sby, str(length) + ' km', transform=tmc,
            horizontalalignment='center', verticalalignment='bottom')
'''
# Given range from a to b
# where a is negative and b is positive
# return an array of tick locations including 0 and spaced interval amounts above and below 0
def calc_ticks(a,b,interval):
    left = np.arange(0, a-1, -interval)
    right = np.arange(0, b+1, interval)
    ticks = np.concatenate([left,right])
    ticks = np.sort(ticks)
    ticks = np.unique(ticks)
    return ticks

# Set the x and y ticks supporting both
# Axis ticks and Grid Lines
def set_ticks(ax_n, interval):
    x0, x1 = ax_n.get_xbound()
    y0, y1 = ax_n.get_ybound()

    yticks, xticks = calc_ticks(y0, y1, interval), calc_ticks(x0, x1, interval)

    ax_n.set_yticks(yticks)
    ax_n.set_xticks(xticks)

    # Not sure we need to set range
    # if some plots look are too small for radar might want to try this...
    #yrange, xrange = (y0, y1), (x0, x1)
    #ax_n.set_ylim(yrange)
    #ax_n.set_xlim(xrange)

# Define how axis ticks are labeled
def set_km_axis_formatter(ax_n):
    # Add km to x-axis.
    def km(x, pos):
        x = int(x/1000)
        return '{} km'.format(x)

    km_formatter = FuncFormatter(km)
    ax_n.xaxis.set_major_formatter(km_formatter)
    ax_n.yaxis.set_major_formatter(km_formatter)
'''   
################################################################################################
##########
# Classes
##########
### set up Master_Plt class (along with subclasses R_Plt (radar plots), and TS_Plt (timeseries plot))
class Master_Plt:
    ##Class Variables
    Domain = 'place holder' # The extent of the area to be plotted
    display = 'place holder' # Define pyart display object for plotting radarfile

    def __init__(self, Data):
        self.Data = Data #(the data dict)

        ## Set some parameters for the plot such as fontsize etc
        #### * * * * * * * * * * * * * * * * * * * * * * * * * *
        #  print(plt.rcParams)
        plt.rc('font', size= 15)         # controls default text sizes
        plt.rc('axes', labelsize= 21) # fontsize of the axes title, and x and y labels
        plt.rc('legend', fontsize= 23, borderpad=.5, facecolor='white', edgecolor= 'black', shadow=True, fancybox=True, framealpha=1)       # legend fontsize
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
                self.ts_gs = GridSpecFromSubplotSpec(len(config.Time_Series), 1, subplot_spec=self.outer_gs[1, :])
            if len(config.Time_Series) == 2:
                self.outer_gs= GridSpec(nrows=2, ncols=1, height_ratios=[4,3], hspace=.1)
                self.ts_gs = GridSpecFromSubplotSpec(len(config.Time_Series), 1, subplot_spec=self.outer_gs[1, :], height_ratios=[2, 3], hspace=0)
            self.r_gs = GridSpecFromSubplotSpec(1, len(config.r_mom), subplot_spec=self.outer_gs[0, :])

        # This is the layout for time series only
        if len(config.Time_Series) != 0 and len(config.r_mom) == 0:
            print('this is the layout for time series only')
            self.fig= plt.figure(figsize=(32,20))

        # This is the layout for radar only
        if len(config.r_mom) != 0 and len(config.Time_Series) == 0:
            print('this is the layout for radar plots only')

        ## Establish some vars that will be helpful to the plot later on
        #### * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
        #if we are plotting radar
        if len(config.r_mom) != 0:
            #####
            #maybe read in the radar data here
            #####
            
            #redefine the classvariable for Domain and display
            Master_Plt.Domain = Platform.getLocation(Data[config.Centered_Pform], offsetkm= config.offsetkm)
            #  self.Domain_Bbox = Bbox.from_extents(self.Domain.xmin, self.Domain.ymin, self.Domain.xmax, self.Domain.ymax)
            
            # Define pyart display object for plotting radarfile
            Master_Plt.display = pyart.graph.RadarMapDisplay(Data['P_Radar'].rfile)
            
            # Set the projection of the radar plot, and transformations 
            self.R_Proj = self.display.grid_projection
            self.trans_Proj = ccrs.PlateCarree()

    # * * * * * * *
    def tick_grid_settings(self, ax, radar=None, ts=None, interval=None, twinax=False):
        ##Axis Limits
        # # # # # # # #
        #if you are calling this for a timeseries subplot set the limits for the plot
        if ts != None:
            # Xaxis 
            # if desired this subsets the timerange that is displayed in the timeseries
            if config.ts_extent != None:
                ax.set_xlim(Platform.Scan_time - timedelta(minutes=config.ts_extent), Platform.Scan_time + timedelta(minutes=config.ts_extent))
            
            # Yaxis 
            if ts in ['Thetav','Thetae']:
                ax.set_ylim(self.Data['Var'].global_min, self.Data['Var'].global_max)
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
            # set up minor ticks (should be a multiple of ten intervals (ie 10,20,30... min spans)
            ax.xaxis.set_minor_locator(AutoMinorLocator(6)) 
            ax.xaxis.set_major_locator(mdates.AutoDateLocator())
            
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

                # Not sure we need to set range
                # if some plots look are too small for radar might want to try this...
                #yrange, xrange = (y0, y1), (x0, x1)
                #ax_n.set_ylim(yrange)
        # + + + + + + + + + + + + ++ + +
        
        ##Tick Display Characteristics
        # # # # # # # # # # # # # # # #
        ax.tick_params(which='minor', axis='both', width=2, length=7, color='grey')
        ax.tick_params(which='major', axis='both', width=2, length=10, color='black')
        # + + + + + + + + + + + + ++ + +

        ##Grid Display Characteristics
        # # # # # # # # # # # # # # # #
        if ts != None:
            ax.tick_params(which='major', axis='both', grid_linewidth=2.5)
            ax.tick_params(which='major', axis='y', grid_color='grey', grid_linewidth=2, grid_linestyle='--', grid_alpha=.8)
            ax.tick_params(which='minor', axis='y', grid_linestyle=':')
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
            #  ax.xaxis.set_visible(False)
        # + + + + + + + + + + + + ++ + +

    # * * * * * * *
    def time_series(self, ts, ax_n, Data, print_long, e_test):
        ''' Plot a time series (sub)plot; could be of pvar info from various intstruments or of wind info
        ----
        INPUTS:
        ts: ........fill in
        Data: dict as described in plotting defn
        print_long & e_test: bool str as described in the ppi defn '''
        if print_long == True: print('~~~~~~~~~~~Made it into time_series~~~~~~~~~~~~~~~~~')

        #  t_TS.T_Plt_settings(ts,ax_n)
        TSleg_elements = [] #empty list to append the legend entries to for each subplot that is actually plotted
        ## MAKE THE TIMESERIES
        #### * * * * * * * * *
        if ts in ['Thetav', 'Thetae']:
            for p in Data.values():
                if isinstance(p, Torus_Insitu):
                    if print_long == True: print('Plotting '+str(p.name)+' on time series')
                    ## If the platform matches a type listed in TS_masking use masked data for the time series; else use the unmasked data
                    if p.type in config.TS_masking[:]: plotting_data= p.ts_mask_df
                    else: plotting_data= p.df[config.p_var].values

                    ## Plot
                    ax_n.plot(p.df['datetime'], plotting_data, linewidth=3, color=p.l_color) #assigning label= is what allows the legend to work
                    TSleg_entry = Line2D([], [], label=p.leg_str, linewidth=12, color=p.l_color)
                    TSleg_elements.append(TSleg_entry)

            ## Set up XY axes tick locations and Labels
            self.tick_grid_settings(ax=ax_n, ts=ts)
            ax_n.set_ylabel(Data['Var'].v_lab)

            ## Plot legend
            leg = ax_n.legend(handles= TSleg_elements, loc='center left')

        # * * *
        if ts == 'Wind':
            p = Data[config.Wind_Pform]
            if print_long == True: print('Plotting '+str(p.name)+' on time series')

            ax_n.plot(p.df['datetime'], p.df['spd'])
            ax_n.fill_between(p.df['datetime'], p.df['spd'], 0)
            ax_n.set_ylabel('Wind Speed')
            self.tick_grid_settings(ax=ax_n, ts=ts)
            TSleg_entry = Line2D([], [], label='Wind Spd', linewidth=12, color='tab:blue')
            TSleg_elements.append(TSleg_entry)

            ax_n.axhline(0, color='k', linewidth=5, zorder=10)

            ax_2 = ax_n.twinx()
            ax_2.plot(p.df['datetime'], p.df['dir'], '.k', linewidth=.05)
            ax_2.set_ylabel('Wind Dir ($^{\circ}$)')
            self.tick_grid_settings(ax=ax_2, ts=ts, twinax=True)
            TSleg_entry = Line2D([], [], marker='.', color='black', label='Wind Dir', markersize=26)
            TSleg_elements.append(TSleg_entry)

            
            ## Plot legend
            leg = ax_n.legend(handles= TSleg_elements, loc='center left')
            #  h, le = ax_n.get_legend_handles_labels()
            #  leg = ax_n.legend(loc='center left')
            leg.set_title(config.Wind_Pform, prop=self.leg_title_font)
            leg.remove()
            ax_2.add_artist(leg)
            
        #if plotting more than one time series then only include the x axis label and ticks to the bottom timeseries
        #  num_of_TS= len(config.Time_Series)
        #  if num_of_TS != 1 and ax_n.rowNum != (num_of_TS-1): plt.setp(ax_n.get_xticklabels(), visible=False)

        ## If makeing the timeseries in conjunction with radar subplots set up vertical lines that indicate the time of
        #  radar scan and the timerange ploted (via filled colorline) on the radar plots
        if len(config.r_mom) != 0:
            ax_n.axvline(Platform.Scan_time, color='r', linewidth=4, alpha=.5)
            ax_n.axvspan(Platform.Scan_time - timedelta(minutes=config.cline_extent), Platform.Scan_time + timedelta(minutes=config.cline_extent), facecolor='0.5', alpha=0.4)

        ##Label Formatter
        # # # # # # # # #
        ax_n.set_xlabel('Time (UTC)')
        ax_n.yaxis.set_label_coords(-.03, .5)

        #  print(ax_n.lines)
        #  print(vars(ax_n))
        if print_long == True: print('~~~~~~~~~~~Made it through time_series~~~~~~~~~~~~~~')


    # * * * * * *  *
    def radar_subplots(self, mom, ax_n, Data, leg, print_long, e_test):
        ''' Plots each of the radar subplots including the marker overlays corresponding to the additional platforms
        ----
        INPUTS:
        mom: str, indicates which radar moment you would like to plot (ie Reflectivity, Velocity, Spectrum Width etc)
        Data: dict as described in the plotting defn
        leg: bool str whether or not you want a legend associated with this particular subplot
        print_long & e_test: bool str as described in the ppi defn
        '''
        if print_long == True: print('~~~~~~~~~~~made it into radar_subplots~~~~~~~~~~~~~~')

        ## SET UP VARS FOR EACH RADAR MOMENTS
        if mom == 'refl':
            p_title, c_label, c_scale = 'Reflectivity', 'Radar Reflectivity [dbz]', 'pyart_HomeyerRainbow'
            if Data['P_Radar'].name in pform_names('KA'):
                field, vminb, vmaxb, sweep = 'refl_fix', -30., 30., Data['P_Radar'].swp
            if Data['P_Radar'].name == 'NOXP':
                field, vminb, vmaxb, sweep = 'DBZ', -30., 30., Data['P_Radar'].swp
            if Data['P_Radar'].name == 'WSR88D':
                field, vminb, vmaxb, sweep =  'reflectivity', -10., 75., Data['P_Radar'].swp[0]
        elif mom == 'vel':
            p_title, c_label, c_scale = 'Radial Velocity', 'Velocity [m/s]', 'pyart_balance'
            if Data['P_Radar'].name in pform_names('KA'):
                field, vminb, vmaxb, sweep = 'vel_fix', -40., 40., Data['P_Radar'].swp
            if Data['P_Radar'].name == 'NOXP':
                field, vminb, vmaxb, sweep = 'VEL', -40., 40., Data['P_Radar'].swp
            if Data['P_Radar'].name == 'WSR88D':
                field, vminb, vmaxb, sweep = 'velocity', -40., 40., Data['P_Radar'].swp[1]
        else: print("Hey what just happened!\n Check the Radar moments for spelling")

        ## Plot the radar
        ax_n.set_title(p_title, y=-.067, fontdict=self.Radar_title_font)
        self.display.plot_ppi_map(field, sweep, ax=ax_n, cmap=c_scale, vmin=vminb, vmax=vmaxb, width=config.offsetkm*2000, height=config.offsetkm*2000, title_flag=False, colorbar_flag=False, embelish=False)
        
        ## PLOT PLATFORMS AS OVERLAYS(ie marker,colorline etc) ON RADAR
        #  iterate over each object contained in dict Data (returns the actual objects not their keys)
        legend_elements = [] #empty list to append the legend entries to for each platfrom that is actually plotted
        for p in Data.values():
            ######
            #Plot Inistu Torus Platforms (if desired and available)
            if isinstance(p, Torus_Insitu):
                if print_long == True: print(p.name)
                legend_elements.append(p.leg_entry)
                self.plot_TORUSpform(p, Data, ax_n, print_long, e_test)

            ######
            #Plot Stationary Inistu Platforms ( aka mesonets and ASOS) if desired and available
            if isinstance(p, Stationary_Insitu):
                if print_long == True: print(p.name)
                # determine if any of the sites fall within the plotting domain
                sites_subdf, valid_sites = p.grab_pform_subset( p, print_long, e_test, Data, bounding= self.Domain)
                # if there are sites within the domain plot the markers and include in the legend
                if valid_sites == True:
                    legend_elements.append(p.leg_entry)
                    ax_n.plot(sites_subdf.lon, sites_subdf.lat, transform=self.trans_Proj, marker=p.m_style, linestyle='None', markersize= p.m_size, color=p.m_color)
                    # include labels for the sites on the plots
                    if p.marker_label == True:
                        for x, y, lab in zip(sites_subdf['lon'], sites_subdf['lat'], sites_subdf['Stn_ID']):
                            trans = self.R_Proj.transform_point(x+.009, y-.002, self.trans_Proj)
                            ax_n.text(trans[0], trans[1], lab, fontsize= 20)

            #####
            #Plot Radar Platforms (if desired and available)
            if isinstance(p, Radar):
                if p.type != 'MAINR':
                    if print_long == True: print(p.name)
                    #det if (any of) the radar is located within the area included in the plot domain at the time of the plot
                    sites_subdf, valid_sites = p.grab_pform_subset(p, print_long, e_test, Data, bounding= self.Domain)

                    #if these conditions are met then plot the radar marker(s)
                    if valid_sites == True:
                        legend_elements.append(p.leg_entry)
                        if p.type == 'KA' or p.type == 'NOXP':
                            ## Plot the marker
                            ax_n.plot(p.lon, p.lat, transform=self.trans_Proj, marker=p.m_style, color=p.m_color, markersize=p.m_size,
                                      markeredgewidth=5, path_effects=[PathEffects.withStroke(linewidth=15, foreground='k')], zorder=10)
                            ## Optional textlabel on plot
                            if p.marker_label == True:
                                ax_n.text(p.lon+.009, p.lat-.002, p.name, transform=trans_Proj, path_effects=[PathEffects.withStroke(linewidth=4, foreground='xkcd:pale blue')])
                            
                            if p.type == 'KA':
                                ## Plot RHI spokes
                                if config.rhi_ring == True: self.rhi_spokes_rings(p)
                            
                        if p.type == 'WSR':
                            ## Plot the marker
                            ax_n.plot(sites_subdf['lon'], sites_subdf['lat'], transform=self.trans_Proj, marker=p.m_style, color=p.m_color, markersize=p.m_size,
                                      markeredgewidth=5, path_effects=[PathEffects.withStroke(linewidth=12, foreground='k')], zorder=8)
                            # include labels for the sites on the plots
                            if p.marker_label == True:
                                for x, y, lab in zip(sites_subdf['lon'], sites_subdf['lat'], sites_subdf['R_Name']):
                                    trans = self.R_Proj.transform_point(x, y, self.trans_Proj)
                                    ax_n.text(trans[0]+.03, trans[1]-.01, lab, fontsize= 27)

        ## PLOT BACKGROUND FEATURES
        if config.country_roads == True:
            ox.config(log_file=True, log_console=True, use_cache=True) #the config in this line has nothing to do with config.py
            G = ox.graph_from_bbox(self.Domain.ymax, self.Domain.ymin, self.Domain.xmax, self.Domain.xmin)
            ox.save_load.save_graph_shapefile(G, filename='tmp'+str(0), folder=config.g_roads_directory , encoding='utf-8')
            fname = config.g_roads_directory + 'tmp'+str(0)+'/edges/edges.shp'
            shape_feature = ShapelyFeature(Reader(fname).geometries(), self.trans_Proj, edgecolor='gray', linewidth=0.5)
            ax_n.add_feature(shape_feature, facecolor='none')
            shutil.rmtree(config.g_roads_directory+'tmp'+str(0)+'/')
        if config.hwys == True:
            fname = config.g_roads_directory+'GPhighways.shp'
            shape_feature = ShapelyFeature(Reader(fname).geometries(), self.trans_Proj, edgecolor='grey')#edgecolor='black')
            ax_n.add_feature(shape_feature, facecolor='none')
        if config.county_lines == True:
            fname = config.g_roads_directory+'cb_2017_us_county_5m.shp'
            shape_feature = ShapelyFeature(Reader(fname).geometries(), self.trans_Proj, edgecolor='gray')
            ax_n.add_feature(shape_feature, facecolor='none', linewidth=1.5, linestyle="--")
        if config.state_lines == True:
            states_provinces = cartopy.feature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='10m', facecolor='none')
            ax_n.add_feature(states_provinces, edgecolor='black', linewidth=2)
        if print_long == True: print('~~~~~~~~~~~Made it through radar_subplots~~~~~~~~~~~')

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
            l = ax_n.legend(handles=legend_elements, loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(-0.1,.5), handlelength=.1)#, title="Platforms")
            #  h, le = ax_n.get_legend_handles_labels()
            #  l = ax_n.legend(loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(0,.5), handlelength=.1)#, title="Platforms")
            l.set_title("Platforms", prop=self.leg_title_font)
        if leg == False:  #this means you are currently making the right subplot
            ## Plot platform colorbar
            #set up colorbar axis that will be as tall and 5% as wide as the 'parent' radar subplot
            cbar_ax = inset_axes(ax_n, width= '5%', height= '100%', loc='center left', bbox_transform=ax_n.transAxes, bbox_to_anchor=(-.29, 0,1,1))
            cbar = plt.colorbar(Data['Var'].CS3, cax=cbar_ax, orientation='vertical', label=Data['Var'].v_lab, ticks=MaxNLocator(integer=True))#,ticks=np.arange(Data['p_var'].global_min, Data['p_var'].global_max+1,2))
        
        ###################
        self.tick_grid_settings(ax=ax_n, radar=True, interval=10*1000)
        scale_bar(ax_n, 10) # 10 KM

        #  set_ticks(ax_n, 10 * 1000) # 10K meters
        #  set_km_axis_formatter(ax_n)

        ax_n.grid(True)
        ###################
        test= ax_n.get_xticks()
        print(test)
        ax_n.set_xticklabels([])    
        #  ax_n.xlocator = MultipleLocator(5000)

        #  ax_n.xaxis.set_major_locator(MultipleLocator(5000)) # set up major tick marks (this is set up to go by 5's will want to change for diff vars)
        #  ax_n.locator_params(MultipleLocator())

        #  print(vars(ax_n))
        #  print(ax_n.lines)
        #  print(ax_n.get_lines)

        if print_long == True: print('~~~~~~~~~~~Made it through radar_subplots~~~~~~~~~~~')

    # * * *
    def plot_TORUSpform(self, pform, Data, ax, print_long, e_test, border_c='xkcd:light grey', labelbias=(0,0)):
        ''' Plot the in situ platform markers, barbs and pathline
        ----
        INPUTS: file: the in situ pandas dataframe
                var: dictionary containing info relating to p_var (ie name, max, min)
                p_attr: dictionary containing platform styling info (color, shape etc)
                ax: axes of the subplot to plot on
        Optional Inputs: border_c, labelbias: color of background for the pathline, if you want to add labels directly to the plot this can offset it from the point
        '''
        if print_long == True: print('made it into platform_plot')

        #grab the subset of data of +- interval around radar scan
        p_sub, p_deploy = pform.grab_pform_subset(pform, print_long, e_test, Data, time_offset=config.cline_extent)

        if p_deploy == False:
            if print_long == True: print('The platform was not deployed at this time')

        #if there is data for the platform that falls within the time and location of interest
        elif p_deploy == True:
            ##Plot the line that extends +/- min from the platform location; The colorfill indicates values of the specifiec p_var (ie Thetae etc)
            #  fill in the values of the colorfill
            C = cmocean.cm.thermal((p_sub[config.p_var].values - Data['Var'].global_min) / (Data['Var'].global_max - Data['Var'].global_min))
            ax = plt.gca()

            #  ax.plot(p_sub.lon, p_sub.lat, c=border_c, linewidth=10.5, transform=ccrs.PlateCarree(), zorder=3)
            #  ax.plot(p_sub.lon, p_sub.lat, c=C, linewidth=7.5, transform=ccrs.PlateCarree(), zorder=4)
            for i in np.arange(len(p_sub['lon']) - 1):
                #This is the border of the colorline
                x, y = p_sub['lon'].values, p_sub['lat'].values
                ax.plot([x[i], x[i+1]], [y[i], y[i+1]], c=border_c, linewidth=10.5, transform=self.trans_Proj, zorder=3)
                #This is the colorramp colorline
                ax.plot([x[i], x[i+1]], [y[i], y[i+1]], c=C[i], linewidth=7.5, transform=self.trans_Proj, zorder=4)
                #  ax.plot([x[i], x[i+1]], [y[i], y[i+1]], c=C[i], linewidth=7.5, transform=ccrs.PlateCarree(), zorder=4,path_effects=[PathEffects.withStroke(linewidth=10.5, color=border_c)] )

            #determine the row of the pform data at scantime
            C_Point=p_sub.loc[p_sub['datetime'].sub(pform.Scan_time).abs().idxmin()]

            #plot the platform marker at the time closest to the scantime (aka the time at the halfway point of the subset platform dataframe)
            ax.plot(C_Point.lon, C_Point.lat, transform=self.trans_Proj, marker=pform.m_style, markersize=pform.m_size,
                    markeredgewidth='3', color=pform.m_color, path_effects=[PathEffects.withStroke(linewidth=12, foreground='k')], zorder=10)
            #plot labels for the marker on the plot itself
            if config.TIn_lab == True:
                plt.text(C_Point.lon+labelbias[0], C_Point.lat+labelbias[1], pform.name, transform=self.trans_Proj, fontsize=20, path_effects=[patheffects.withstroke(linewidth=4)])

            #plot a dot at the end of the colorline in the direction the platform is  moving (aka the last time in the subset dataframe)
            ax.plot(p_sub['lon'].iloc[-1], p_sub['lat'].iloc[-1], transform=self.trans_Proj, marker='.', markersize=10, markeredgewidth='3', color='k')

            #plot windbarbs
            #  p_sub.iloc[::x,col_index] returns every x'th value
            stationplot = StationPlot(ax, p_sub.loc[::30, 'lon'], p_sub.loc[::30, 'lat'], transform=self.trans_Proj)
            stationplot.plot_barb(p_sub.loc[::30, 'U'], p_sub.loc[::30, 'V'], sizes=dict(emptybarb=0), length=7)

        if print_long == True: print('made it through platform_plot')

    # * * *
    def rhi_spokes_rings(self, pform):
        ''' Plot the RHI spoke and ring for a radar
        '''
        if print_long == True: print('made it into rhi_spokes_rings')

        # if there is not actually rhi info then it will not plot a ring and not stop the code
        if np.isnan(pform.rhib) == True or np.isnan(pform.rhie)== True: print('Could not plot RHI spokes')

        #produce spoke and ring
        else:
            for j in range(int(pform.rhib), int(pform.rhie)+1, 10):
                ang = pform.head + j
                if ang > 360.: ang= int(ang - 360.)
                #  radius = Data['P_Radar'].rfile.range['data'][-1]-500.)/1000.
                if pform.type == 'KA': radius= 20.905
                else: print('code not written for other radars yet')

                #this plots a circle that connects the spokes
                latArray, lonArray = [], []
                for bearing in range(int(pform.head + pform.rhib), int(pform.head + pform.rhie+1)): #degrees of sector
                    lat2, lon2 = Platform.getLocation(pform, radius, given_bearing=bearing)
                    latArray.append(lat2)
                    lonArray.append(lon2)
                self.display.plot_line_geo(lonArray, latArray, marker=None, color='grey', linewidth=.25) #this plots a circle that connects the spokes

                #plt the spokes
                C, D = Platform.getLocation(pform, radius, given_bearing = ang)
                self.display.plot_line_geo([pform.lon, D], [pform.lat, C], marker=None, color='k', linewidth=0.5, linestyle=":")

                ## optional labels
                if config.RHI_lab == True:
                    if np.logical_and(C>ymin, np.logical_and(C<ymax, np.logical_and(D>xmin, D<xmax))):
                        plt.text(D, C, str(ang), horizontalalignment='center', transform=self.trans_Proj, zorder=9,
                                 path_effects=([PathEffects.withStroke(linewidth=4, foreground='xkcd:pale blue')]))
            if print_long == True: print('made it through rhi_spokes_rings')

#######################################################
        #  plt.rc('axes', xmargin = 0,  ymargin=0)
        #  if ts == 'Wind':
            #  with mpl.rc_context(rc={})

    #  def things_for_rc():
        #  plt.rc('font',weight='bold')
        #  plt.rc('xtick.major', size=5, pad=7)
        #  plt.rc('xtick', labelsize=15)
        #  plt.rc('grid',c='.5', ls='-',lw=5)
#
    #  @ticker.FuncFormatter
    #  def ticklab_format(x, pos):
        #  return f'[{x:.2f}]'
    #  ax_n.xaxis.set_major_formatter(ticklab_format)

    #  @staticmethod
    #  def plot_bground_features():
##########

##########################################################################
###########################
## PLOTTING DEFINTIONS  ###
###########################
def plotting(Data, print_long, e_test, start_comptime):
    ''' Initial plotting defenition: sets up fig size, layout, font size etc and will call timeseries and radar subplots
    ----------
    INPUTS:
    Data: (dictionary)
        contains objects corresponding to all available datasets (radar, torus insitue platforms etc) and  which particular variable
        should be plotted on time series and colorlines
    print_long & e_test: (bool strings)
        True/False vars that control how much information is printed out to the terminal window. Helpful for debugging but can be overkill
    start_comptime:
        Time at which you first begain plotting this particular image (will help to report out how long it took to create image)
    '''
    if print_long == True: print('~~~~~~~~~~~Made it into plotting~~~~~~~~~~~~~~~~~~~~~')
    #initilize the plot object(will have info about plotlayout and such now)
    PLT = Master_Plt(Data)

    ## Make the Radar Subplots (if applicable)
    #### * * * * * * * * * * * * * * * * * ** * * * *
    if len(config.r_mom) != 0:
        for (subcol,), mom in np.ndenumerate(config.r_mom):
            print('Radar plot: '+ mom +', Outer GridSpec Pos: [0, :], SubGridSpec Pos: [:, '+str(subcol)+']')
            ax_n= PLT.fig.add_subplot(PLT.r_gs[:, subcol], projection= PLT.R_Proj)
            if subcol==0: leg = True
            else: leg = False
            PLT.radar_subplots(mom, ax_n, Data, leg, print_long, e_test)

    ## Make the Times Series Subplots (if applicable)
    #### * * * * * * * * * * * * * * * * * ** * * * *
    if len(config.Time_Series) != 0:
        for (subrow,), ts in np.ndenumerate(config.Time_Series):
            print('Time Series plot: '+ ts +', Outer GridSpec Pos: [1, :], SubGridSpec Pos: ['+str(subrow)+', :]')
            if subrow==0: ax_n= PLT.fig.add_subplot(PLT.ts_gs[subrow,:])
            else:  ax_n=PLT.fig.add_subplot(PLT.ts_gs[subrow,:], sharex=ax_n)
            PLT.time_series(ts, ax_n, Data, print_long, e_test)

    ## Finish plot
    #### * * * * *
    #Add a title in the plot
    if Data['P_Radar'].name in pform_names('KA') or Data['P_Radar'].name =='WSR88D':
        plt.suptitle(Data['P_Radar'].site_name+' '+str(config.p_tilt)+r'$^{\circ}$ PPI '+Data['P_Radar'].fancy_date_str, y=.92)
    else: 
        plt.suptitle(Data['P_Radar'].name+' '+str(config.p_tilt)+r'$^{\circ}$ PPI '+Data['P_Radar'].fancy_date_str, y=.92)
    file_string = '_'.join(config.Time_Series)
    
    #This is the directory path for the output file
    outdir_name = config.g_plots_directory+config.day+'/plots/'+Data['P_Radar'].dir_name+'/'+config.p_var+'/'+str(config.p_tilt)+'_deg/'
    # Setup Function Cache for speedup, if path directory does not make the directory 
    if not os.path.exists(outdir_name): Path(outdir_name).mkdir(parents=True)
    
    #add the file name 
    if Data['P_Radar'].name in pform_names('KA'):
        output_name= outdir_name+Data['P_Radar'].site_name+'_'+Platform.Scan_time.strftime('%m%d_%H%M')+'_'+file_string+'.png'
    if Data['P_Radar'].name == 'NOXP':
        output_name = outdir_name+Platform.Scan_time.strftime('%m%d_%H%M')+'_'+file_string+'_NOXP.png'
    if Data['P_Radar'].name == 'WSR88D':
        output_name = outdir_name+Platform.Scan_time.strftime('%m%d_%H%M')+'_'+file_string+'_'+Data['P_Radar'].site_name+'.png'
    print(output_name)

    plt.savefig(output_name, bbox_inches='tight', pad_inches=.3)
    timer(start_comptime, time.time())
    plt.close()

    print('\a') ## Makes a ding noise
    if print_long == True: print('~~~~~~~~~~~made it through plotting~~~~~~~~~~~~~~~~~~')
    print('Done Plotting \n \n***************************************************************************************************')

##############################################################################################
#######################
## RADAR DEFINITIONS ##
#######################
# * * *
def det_nearest_WSR(p_df):
    ''' locate the nearest WSR88D site to the specified insitu instruments
    '''
    #find the locations of all WSR88D sites (outputs dict in format {Site_ID:{lat:...,lon:...,elav:...], ...})
    all_WSR = pyart.io.nexrad_common.NEXRAD_LOCATIONS
    #  print(json.dumps(all_WSR_locs, sort_keys=True, indent=4))

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
    return r_ofintrest

# * * *
def get_WSR_from_AWS(start, end, radar_id, download_directory):
    ''' Retrieve the NEXRAD files that fall within a timerange for a specified radar site from the AWS server
    ----------
    INPUTS
    radar_id : string,  four letter radar designation
    start: datetime,  start of the desired timerange
    end: datetime,  end of the desired timerange
    download_directory: string,  location for the downloaded radarfiles
    -------
    RETURN
    radar_list : Py-ART Radar Objects
    '''
    # Create this at the point of use # Otherwise it saves everything and eventually crashes
    conn = nexradaws.NexradAwsInterface()

    #Determine the radar scans that fall within the time range for a given radar site
    scans = conn.get_avail_scans_in_range(start, end, radar_id)
    print("There are {} scans available between {} and {}\n".format(len(scans), start, end))

    #download the files that were identified
    #results = conn.download(scans[0:4], filesys+'TORUS_Data/'+day+'/radar/Nexrad/Nexrad_files/', keep_aws_folders=False)
    #results = conn.download(scans, filesys+'TORUS_Data/'+day+'/radar/Nexrad/Nexrad_files/', keep_aws_folders=False)
    #results = conn.download(scans[0:4], temploc+day+'/radar/Nexrad/Nexrad_files/', keep_aws_folders=False)

    # Don't download files that you already have...
    path =  download_directory + config.day +'/radar/Nexrad/Nexrad_files/'

    if not os.path.exists(path): Path(path).mkdir(parents=True)

    # missing_scans is a list of scans we don't have and need to download
    # create_filepath returns tuple of (directory, directory+filename)
    # [-1] returns the directory+filename
    missing_scans = list(filter(lambda x: not Path(x.create_filepath(path,False)[-1]).exists(), scans))

    # missing files is the list of filenames of files we need to down load
    missing_files = list(map(lambda x: x.create_filepath(path,False)[-1], missing_scans))
    print("missing ", len(missing_files), "of ", len(scans), " files")
    print(missing_files)

    results = conn.download(missing_scans, path, keep_aws_folders=False)
    print(results.success)
    print("{} downloads failed: {}\n".format(results.failed_count,results.failed))
    #print("Results.iter_success : {}\n".format(results.iter_success()))

    # missing_scans_after is a list of scans we don't have (download failed)
    # create_filepath returns tuple of (directory, directory+filename)
    # [-1] returns the directory+filename
    missing_files_after = list(filter(lambda x: not Path(x.create_filepath(path,False)[-1]).exists(), scans))

    if len(missing_files_after) > 0:
        print("ERROR: Some Radar Scans Missing")
        print(missing_files_after)
        exit()

    radar_files = list(map(lambda x: x.create_filepath(path,False)[-1], scans))

    # Return list of files
    return radar_files

# * * *
def read_from_nexrad_file(radar_file):
    radar = pyart.io.read_nexrad_archive(radar_file)
    return radar
def read_from_KA_file(radar_file):
    radar = pyart.io.read(radar_file)
    return radar
def read_from_NOXP_file(radar_file):
    #  radar = pyart.io.read(radar_file)
    radar = pyart.io.read_cfradial(radar_file)
    return radar
# Note: Cached version is cached on the file name, not the file contents.
# If file contents change you need to invalidate the cache or pass in the file contents directly to this function
# function_cache_memory = Memory(config.temploc, verbose=1)
function_cache_memory = Memory(config.g_cache_directory,verbose=1)
cached_read_from_nexrad_file = function_cache_memory.cache( read_from_nexrad_file )
cached_read_from_KA_file = function_cache_memory.cache( read_from_KA_file )
cached_read_from_NOXP_file = function_cache_memory.cache( read_from_NOXP_file )


# * * *
def plot_radar_file(r_file, Data, subset_pnames, print_long, e_test, swp_id= None):
    print("open_pyart, scan file_name = {}\n".format(r_file))
    start_comptime = time.time()

    if config.Radar_Plot_Type == 'WSR_Plotting':
        #open file using pyart
        try: radar = cached_read_from_nexrad_file(r_file)
        except: 
            radar = 'NO FILE' 
            print("Failed to convert file: "+str(r_file))
        valid_time = True
    
    if config.Radar_Plot_Type == 'KA_Plotting':
        head_tail= os.path.split(r_file)
        time_string= str(20)+str(head_tail[1][13:24])
        rtime = datetime.strptime(time_string, "%Y%m%d%H%M%S")
        print(rtime)
        valid_time = time_in_range(config.tstart, config.tend, rtime)
        print(valid_time)
        #  in_time_range= True
        #  if config.tstart != None:
            #  if rtime >= config.tstart: pass
            #  else: in_time_range = False
        #  if config.tend != None:
            #  if rtime <= config.tend: pass
            #  else: in_time_range = False
        
        #  if in_time_range == False: print(rfile+' was not in the timerange being plotted')
        if valid_time == False: print(r_file+' was not in the timerange being plotted')
        else:
            ## Read the radar file
            radar = cached_read_from_KA_file(r_file)
            ## If file contains a ppi scan proceed; we are currently not interested in plotting the RHI scans
            if radar.scan_type == 'ppi':
                for swp_id in range(radar.nsweeps):
                    tilt_ang = radar.get_elevation(swp_id) ## Det the actual tilt angle of a given sweep (returns an array)
                    ## Check to see if the radarfile matches the elevation tilt we are interested in
                    if np.around(tilt_ang[0], decimals=1) == config.p_tilt:
                        print("\nProducing Radar Plot:")
                        #  Assign radar fields and masking
                    #  det_radar_fields(radar)
    
    if config.Radar_Plot_Type == 'NOXP_Plotting':
        ## Read the radar file
        radar = cached_read_from_NOXP_file(r_file)
        for swp_id in range(radar.nsweeps):
            tilt_ang = radar.get_elevation(swp_id) ## Det the actual tilt angle of a given sweep (returns an array)
            ## Check to see if the radarfile matches the elevation tilt we are interested in
            if np.around(tilt_ang[0], decimals=1) == config.p_tilt:
                print("\nProducing Radar Plot:")
                #  Assign radar fields and masking
        valid_time = True

    #  try:
    if config.print_radar_info== True: print(radar.info(level='compact'))

    ## Read in radar data and add to Data dict
    ##### + + + + + + + + + + + + + + + + + + +
    if radar != 'NO FILE':
        #  Establish info for the main plotting radar (aka scantime etc) & and locations for other radars (if deployed)
        Data, subset_pnames = Add_to_DATA('RADAR', Data, subset_pnames, print_long, MR_file=radar, swp=swp_id)
        if print_long == True: print(str(Data)+'\n')

        ## Proceed to plot the radar
        ##### + + + + + + + + + + + +
        plotting(Data, print_long, e_test, start_comptime)
        #  except:
            #  print('something went wrong')
    '''if config.Radar_Plot_Type == 'KA_Plotting' and valid_time == False: pass
    else:
        if config.print_radar_info== True: print(radar.info(level='compact'))

        ## Read in radar data and add to Data dict
        ##### + + + + + + + + + + + + + + + + + + +
        #  Establish info for the main plotting radar (aka scantime etc) & and locations for other radars (if deployed)
        Data, subset_pnames = Add_to_DATA('RADAR', Data, subset_pnames, print_long, MR_file=radar, swp=swp_id)
        if print_long == True: print(str(Data)+'\n')

        ## Proceed to plot the radar
        ##### + + + + + + + + + + + +
        plotting(Data, print_long, e_test, start_comptime)
    '''
##############################################################################################

#  parser = argparse.ArgumentParser(description='Process Radar')
# parser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')
#  parser.add_argument("--day", help='override day in config file')

#  args = parser.parse_args()

#  print("args: ", args)
#  if args.day: config.day = args.day
#  print ("using day: ", config.day)

##############################################################################################
#####################################################
# Read in Data that does not change for each image ##
#####################################################
print('\nRead in Torus Platforms')
#print(pform_names('ALL')) #list of ALL possible platforms
subset_pnames = [] #array of platforms that actually have data for the time range we are plotting
Data = {} #dictionary in which all the class objects will be stored (will contain platform data, locations, plotting variables etc)

## Read in the data for the TORUS Insitu platforms (if available)
Data, subset_pnames = Add_to_DATA('TInsitu', Data, subset_pnames, print_long)
## Establish info for the plotting variable (aka max min etc)
Data, subset_pnames = Add_to_DATA('PVar', Data, subset_pnames, print_long)

print('\nRead in Stationary Platforms Arrays')
## Read in the data for the Stationary Array platforms (if available)
Data, subset_pnames = Add_to_DATA('STN_I', Data, subset_pnames, print_long)


###################################
# Create Plot for each radarfile ##
###################################
## If Radar will be plotted
if config.r_plotting == True:
    print('\nYes Plot Radar \n')

    # * * *
    if config.Radar_Plot_Type == 'KA_Plotting':
        ## Get radar files
        path = config.g_mesonet_directory + config.day+'/radar/TTUKa/netcdf/*/dealiased_*'
        radar_files = sorted(glob.glob(path))
        if len(radar_files)==0:
            path = '/Users/severe2/Research/TORUS_Data/' +config.day+'/radar/TTUKa/netcdf/*/dealiased_*'
            radar_files = sorted(glob.glob(path))
            if len(radar_files)==0: print(path)
        ## Proceed to plot the radar
        ##### + + + + + + + + + + + +
        Parallel(n_jobs=config.nCPU, verbose=10)(delayed(plot_radar_file)(r_file, Data, subset_pnames, print_long, e_test) for r_file in radar_files)
    
    # * * * 
    if config.Radar_Plot_Type == 'NOXP_Plotting':
        ## Get radar files
        path = config.g_mesonet_directory + config.day+'/radar/NOXP/'+config.day+'/*/sec/*'
        radar_files = sorted(glob.glob(path))
        
        ## Proceed to plot the radar
        ##### + + + + + + + + + + + +
        Parallel(n_jobs=config.nCPU, verbose=10)(delayed(plot_radar_file)(r_file, Data, subset_pnames, print_long, e_test) for r_file in radar_files)

    # * * *
    if config.Radar_Plot_Type == 'WSR_Plotting':
        #Det the unique radar sites to be plotted
        unique_r_sites=det_nearest_WSR( Data[config.Centered_Pform].df)
        if print_long == True: print(unique_r_sites)

        #set up empty dataframe
        tranges_each_r = pd.DataFrame()
        for Rad_site in unique_r_sites:
            print(Rad_site)
            #  print(Data[config.Centered_Pform].df.iloc[-1])
            #  print(Data[config.Centered_Pform].df[-1])
            trange_r = Data[config.Centered_Pform].df.loc[Data[config.Centered_Pform].df.Radar_ID == Rad_site, ['datetime']].rename(columns={'datetime': Rad_site})
            trange_r_start, trange_r_end = trange_r.min(), trange_r.max()
            tranges_each_r = pd.concat([tranges_each_r, trange_r], axis=1)
            #  print(tranges_each_r)

            print("start "+str(trange_r_start[Rad_site])+ "\nend "+str(trange_r_end[Rad_site])+ "\n ***")
            radar_files = get_WSR_from_AWS(trange_r_start[Rad_site], trange_r_end[Rad_site], Rad_site, config.g_download_directory)
            print('********\n Radar files to process:\n'+ str(radar_files))

            #Hard code the swp numbers that will be associated with a given tilt angle
            if config.p_tilt == .5: swp_id=[0 , 1]
            elif config.p_tilt == 1: swp_id=[2 , 3]
            elif config.p_tilt == 1.5: swp_id=[4 , 5]
            else: print('The tilt angle {} is not hard coded yet for WSR'.format(config.p_tilt))

            #open the downloaded files as pyart objects

            ## Proceed to plot the radar
            ##### + + + + + + + + + + + +
            Parallel(n_jobs=config.nCPU, verbose=10)(delayed(plot_radar_file)(r_file, Data, subset_pnames, print_long, e_test, swp_id= swp_id) for r_file in radar_files)

        print(tranges_each_r)


################################
# Create Timeseries only plot ##
################################
#Only plot timeseries (this code isn't fully fleshed out but in theroy this code is built in such a way to allow for this)
if config.r_plotting == False and config.t_plotting == True:
    print("Plot Timeseries only \n"+ str(config.g_mesonet_directory+config.day+'/mesonets/NSSL/*.nc'))
    time_series(Data)
    fig.savefig('test2.png')
    plt.close()

###########################
plt.rcdefaults()
timer(totalcompT_start, time.time(), total_runtime=True)
print("ALL FINISHED")
