
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
from datetime import datetime, timedelta
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
import itertools
from netCDF4 import num2date
import logging
import csv 
import os
import os.path
import time
import datetime as dt
import scipy
from pathlib import Path
from joblib import Memory, Parallel, delayed
from scipy import ndimage, interpolate
import cmocean
import gc
import glob
import nexradaws
import pyart
import gc
from pympler.asizeof import asizeof
from collections import namedtuple
#this is the file with the plotting controls to access any of the vars in that file use config.var
import config as plot_config
import copy
import singledop
import pynimbus as pyn 
import tracemalloc
#import shapefile
import math 
from shapely.geometry import LineString, Point, Polygon
import seaborn as sns
#appnope is here to prevent slow down of pop up gui
import appnope
appnope.nope()
#  sns.set()
tracemalloc.start()
print(matplotlib.get_backend())
import matplotlib.style as mplstyle
mplstyle.use('fast')
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
from radar_defns import read_from_radar_file, radar_fields_prep, det_nearest_WSR, sweep_index, get_WSR_from_AWS, info_radar_files
from surge_defns import clicker_defn, make_csv, csv_reader, surge_radarbins, meso_radarbins

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

        if config.overlays['Colorline']['Pert'] == True: 
            lab, GMin, GMax = lab+' Pert', -10, 10
            CSCALE= 'RdBu_r'
        else: 
            CSCALE= 'gnuplot'
        self.Tvar_lab, self.Tvar_GMin, self.Tvar_GMax = lab, GMin, GMax

        #settings for the thermo var being plotted on radar subplots
        ## Make a dummy plot to allow for ploting of colorbar
        self.norm =plt.Normalize(GMin, GMax)

        Z = [[0,0],[0,0]]
        levels = np.arange(GMin, GMax+1, 1)
        self.CS3 = plt.contourf(Z, levels, cmap=CSCALE)
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
    def __init__(self, day, Data, config, tilt=None):
        self.Data = Data #(the data dict)
        self.config = config

        ## Set some parameters for the plot such as fontsize etc
        #### * * * * * * * * * * * * * * * * * * * * * * * * * *
        #  print(plt.rcParams)
        plt.rc('font', size= 13)      # controls default text sizes
        plt.rc('axes', labelsize= 12) # fontsize of the axes title, and x and y labels
        plt.rc('xtick', labelsize= 10) # fontsize of the xtick labels
        plt.rc('ytick', labelsize= 10) # fontsize of the ytick labels
        plt.rc('legend', fontsize= 13, borderpad=.45, facecolor='white', 
               edgecolor= 'black', shadow=True, fancybox=True, framealpha=1)# legend fontsize
        plt.rc('figure', titlesize= 20, titleweight='semibold',facecolor='white')  # fontsize of the figure title
        self.leg_title_font={'size':15, 'weight':'bold'}
        self.Radar_title_font= {'fontsize':13}
        self.Contour_label_font= {'fontsize':10}


        ## Establish the Plot size and layout
        #### * * * * * * * * * * * * * * * * *
        def radar_layout(config):
            if 'zoom' in config.r_mom: 
                Num_of_rplts= len(config.r_mom)+self.num_of_surges-1
            else:
                Num_of_rplts= len(config.r_mom)

            #if you specify a start or end time it will be assigned here otherwise will be set to none (full data set)
            if config.radar_controls['layout']['Rows'] == None:
                R_rows = math.trunc(Num_of_rplts / 2)
            if config.radar_controls['layout']['Cols'] == None:
                if ((Num_of_rplts / 2) % 2) in [0,1]: R_cols = 2
                else: R_cols = 3

            if config.radar_controls['layout']['Rows'] != None:
                R_rows = config.radar_controls['layout']['Rows']
                R_cols = round(Num_of_rplts / R_rows) 

            if config.radar_controls['layout']['Cols'] != None:
                R_cols = config.radar_controls['layout']['Cols']
                R_rows = round(Num_of_rplts / R_cols) 
            
            return R_rows, R_cols

        # * * * 
        def lineplot_layout(config):
            if 'clicker' in config.Line_Plot: 
                if self.num_of_surges ==0:
                    L_Rows= len(config.Line_Plot)
                else:
                    L_Rows= len(config.Line_Plot)+self.num_of_surges-1
                L_wratio = [.25,2,2,1]
                L_Cols=4
            else:
                L_Rows= len(config.Line_Plot)
                if config.lineplt_control['Deriv']== True: L_Cols, L_wratio = 2, [2,1]
                else:  L_Cols, L_wratio = 1, [1]
            
            if L_Rows == 1:
                L_hratio=[1] #the relative height alocated to each timeseries [within the space allocated for all timeseries]
                h_space=0
            elif L_Rows == 2:
                if 'Wind' in config.Line_Plot: L_hratio=[2,3]
                else: L_hratio=[1,1]

                if 'clicker' in config.Line_Plot: h_space=.2
                else: h_space= 0

            return L_Rows, L_Cols, L_hratio, L_wratio, h_space

        # * * * 
        def outer_layout(config): 
            if len(config.Line_Plot) == 1:
                outer_hratio = [4,2] #the height ratio between the space alocated for radar imagaes and timeseries
            if len(config.Line_Plot) == 2 or self.num_of_surges > 1:
                outer_hratio = [4,3] 
            return outer_hratio

        # * * * 
        def figsize_options(config):
            #, r_rows):
            self.title_spacer=.93
            #  xfigsize= 15
            #  xfigsize= 20
            xfigsize=self.R_cols*7
            yfigsize= self.R_rows*12#9
            if self.L_Rows > 1:
                yfigsize = yfigsize+ ((self.L_Rows-1)*3)
            self.fig= plt.figure(figsize=(xfigsize,yfigsize))
            

        #######
        if config.Surge_controls['Feature_IDing']['Existing_pnts']['Read_in_surges'] == True: 
            self.num_of_surges, self.surge_ids, self.valid_surges = csv_reader(config, Data, day, tilt, 'Surge')
            if self.num_of_surges !=0:
                self.Surges, self.gate_x, self.gate_y, self.gate_lon, self.gate_lat = surge_radarbins(self, config, Data, Data['P_Radar'].swp)  
        else: 
            self.num_of_surges=0

        if config.Surge_controls['Feature_IDing']['Existing_pnts']['Read_in_mesos'] == True: 
            self.num_of_mesos, self.meso_ids, self.valid_mesos = csv_reader(config, Data, day, tilt, 'Meso')
            if self.num_of_mesos!=0:
                self.Mesos=meso_radarbins(self, day, config, Data, Data['P_Radar'].swp)
        else: 
            self.num_of_mesos=0
            self.meso_ids=None

        #######
        #both ts and r
        if len(config.Line_Plot) != 0 and len(config.r_mom) != 0:
            R_and_L_plt=True
            outer_hratio=outer_layout(config)
            self.outer_gs= GridSpec(nrows=2, ncols=1, height_ratios=outer_hratio, hspace=.17)
        else:
            R_and_L_plt=False

           #######

        if len(config.Line_Plot) == 0:
            self.L_Rows, self.L_Cols= 0, 0 
        elif len(config.Line_Plot) != 0:
            self.L_Rows, self.L_Cols, L_hratio, L_wratio, h_space = lineplot_layout(config)
            if R_and_L_plt == True:
                self.l_gs = GridSpecFromSubplotSpec(self.L_Rows, self.L_Cols, subplot_spec=self.outer_gs[1, :], height_ratios=L_hratio, width_ratios=L_wratio, hspace=h_space)
            else:
                self.l_gs = GridSpec(self.L_Rows, self.L_Cols, height_ratios=L_hratio, width_ratios=L_wratio, hspace=h_space)

        
        #if we are plotting radar
        if len(config.r_mom) == 0:
            self.R_rows, self.R_cols = 0, 0
        elif len(config.r_mom) != 0:
            self.R_rows, self.R_cols = radar_layout(config)
            if R_and_L_plt == True:
                self.r_gs = GridSpecFromSubplotSpec(self.R_rows, self.R_cols, subplot_spec=self.outer_gs[0, :], wspace=0)
            else:
                self.r_gs = GridSpec(self.R_rows, self.R_cols, hspace=0)
            
            ## Establish some vars that will be helpful to the plot later on
            #### * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
            if len(config.r_mom) != 0:
                # Define pyart display object for plotting radarfile
                self.display = pyart.graph.RadarMapDisplay(Data['P_Radar'].rfile)

                ### The extent of the area to be plotted
                self.Domain = Platform.getLocation(Data[config.radar_controls['Centered_Pform']], scan_time=Data['P_Radar'].Scan_time, 
                                                   offsetkm= config.radar_controls['offsetkm'][config.Radar_Plot_Type])
                #  self.Domain_Bbox = Bbox.from_extents(self.Domain.xmin, self.Domain.ymin, self.Domain.xmax, self.Domain.ymax)

                # Set the projection of the radar plot, and transformations
                self.R_Proj = self.display.grid_projection
                self.Proj = ccrs.PlateCarree()

        figsize_options(config)

    # * * * * * * *
    # * * * * * * *
    def tick_grid_settings(self, ax, scan_time, radar=None, ts=None, interval=None, twinax=False):
        ##Axis Limits
        # # # # # # # #
        #if you are calling this for a timeseries subplot set the limits for the plot
        if ts != None:
            # Xaxis
            # if desired this subsets the timerange that is displayed in the timeseries
            if scan_time != None and self.config.lineplt_control['Axis_Extents']['X'] != None:
                ax.set_xlim(scan_time - timedelta(minutes=self.config.lineplt_control['Axis_Extents']['X']),
                            scan_time + timedelta(minutes=self.config.lineplt_control['Axis_Extents']['X']))

            # Yaxis
            if ts in ['Thetav','Thetae']:
                if self.config.lineplt_control['Axis_Extents']['Y'] != None: ax.set_ylim(335,350)
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
            #  plt.setp(ax.get_yticklabels(), labellocation='Center')

        # + + + + + + + + + + + + ++ + +
        ##Tick Display Characteristics
        # # # # # # # # # # # # # # # #
        ax.tick_params(which='minor', axis='both', direction='in', width=2, length=7, color='grey')
        ax.tick_params(which='major', axis='both', direction='in', width=2, length=10, color='black')
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
    
    # * * * * * * * * * * * * * * * * * * * * * * * *
    # * * * * * * * * * * * * * * * * * * * * * * * *
    def radar_subplots(self, P_Radar, mom, day, ax_n, fig, TVARS, Data, leg, thermo_cbar, Surge_ID, Meso_ID, ZOOM=False):
        ''' Plots each of the radar subplots including the marker overlays corresponding to the additional platforms
        ----
        INPUTS: mom: str, indicates which radar moment you would like to plot (ie Reflectivity, Velocity, Spectrum Width etc)
                Data: dict as described in the plotting defn
                leg: bool str whether or not you want a legend associated with this particular subplot
        '''
        if self.config.print_long == True: print('~~~~~~~~~~~made it into radar_subplots~~~~~~~~~~~~~~')

        if isinstance(P_Radar, Radar):
            sweep, rfile, Rname, Scan_time, DISPLAY = P_Radar.swp, P_Radar.rfile, P_Radar.name, P_Radar.Scan_time, self.display
        else:
            sweep, rfile, Rname, Scan_time = P_Radar['swp'], P_Radar['rfile'], P_Radar['name'], P_Radar['Scan_time']
            DISPLAY= pyart.graph.RadarMapDisplay(rfile)

        if Rname == 'WSR88D':
           if mom in ['refl','refl_unmasked', 'refl_despeck']:
               sweep = sweep[0]
           else: 
               sweep = sweep[1]

        ## SET UP VARS FOR EACH RADAR MOMENTS
        #####################################
        testing_testing= False
        testing_bound=None
        #  if mom == 'refl':
            #  mom = 'Ray_angle'
        if mom in ['refl','refl_unmasked', 'refl_despeck']:
            c_scale = 'pyart_HomeyerRainbow'
            if Rname == 'WSR88D':
                vminb, vmaxb = -10., 75.
            else:
                if Rname in pform_names('KA'): vminb, vmaxb= -30., 30.
                elif Rname == 'NOXP': vminb, vmaxb= -10, 60
            #****
            if mom == 'refl':
                c_label= 'Reflectivity [dbz]' 
                if Rname == 'WSR88D': field = 'reflectivity'
                else: field= 'refl_fix'
            if mom == 'refl_unmasked':
                c_label= 'Reflectivity [dbz]' 
                if Rname == 'WSR88D': field = 'reflectivity'
                elif Rname in pform_names('KA'): field= 'reflectivity'
                elif Rname =='NOXP': field= 'DBZ'
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
            vminb,vmaxb= -.15, .15  
        
        ###
        elif mom == 'zoom':
            c_scale = cmocean.cm.balance
            if Rname == 'WSR88D':
               vminb, vmaxb= -40., 40.
            else:
               vminb, vmaxb= -30., 30.
               #  vminb, vmaxb, sweep = -30., 30., Data['P_Radar'].swp
            
            c_label= 'Radial Velocity Smoothed [m/s]'
            field= 'vel_smooth'

        
        ###
        elif mom in ['vel', 'vel_unmasked', 'vel_unfixed', 'vel_despeck', 'vel_smooth', 'vel_savgol']: 
            c_scale = cmocean.cm.balance
            if Rname == 'WSR88D':
               vminb, vmaxb= -40., 40.
            else:
               vminb, vmaxb= -30., 30.
               #  vminb, vmaxb, sweep = -30., 30., Data['P_Radar'].swp
            #****
            if mom == 'vel':
                c_label = 'Radial Velocity [m/s]'
                if Rname == 'WSR88D': field= 'velocity'
                else: field= 'vel_fix'
            elif mom == 'vel_unmasked':
                c_label = 'Radial Velocity [m/s]'
                if Rname == 'WSR88D': field= 'velocity'
                elif Rname in pform_names('KA'): field= 'corrected_velocity'
                elif Rname =='NOXP': field= 'corrected_velocity'
            elif mom == 'vel_unfixed':
                c_label= 'Velocity Unfixed [m/s]'
                if Rname in pform_names('KA'): field= 'velocity'
                elif Rname == 'NOXP': field = 'VEL' 
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
            field, vminb, vmaxb= 'specw_fix', 0., 4
        ###
        elif mom == 'Meso_azi':
            c_label, c_scale = 'Azimuthal Shear 'r' (s$^{-1}$)', 'RdBu'
            field, vminb, vmaxb= 'Meso_azi', -10, 10 
        ###
        elif mom in ['vel_texture','vel_texture_smoothed', 'vel_texture_dea', 'vel_test']:
            c_scale = 'bone_r'#'PuBuGn'#'BuPu'#cmocean.cm.dense#cmocean.cm.speed#Greys #'pyart_balance'
            #  c_scale = 'Greys'#'pyart_balance'
            vminb, vmaxb=  0., 3
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
                #  testing_bound=[0,.5,1,1.5,2,2.5,3,3.5]
                testing_bound=[0,1,2,3,4]
                #  testing_bound=[0,.5,1,1.5,2,2.5]
                vminb, vmaxb=  0., 3
                #  testing_bound=[0,1,1.25, 1.5, 1.75,2,3,]
                c_label= 'Radial Velocity masking test'
                field= 'vel_test'

        ###
        elif mom in ['vel_grad', 'vel_savgol_axis1', 'vel_savgol_axis0']:
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
            if mom == 'az_shear':
                c_label, c_scale = 'Azimulathal Shear (s$^{-1}$)', 'RdBu'
                field, vminb, vmaxb = 'az_shear', -15, 15 
        ###
        elif mom == 'diff':
            p_title, c_label, c_scale = 'Difference (fixed - original)', 'Difference', cmocean.cm.balance#'pyart_balance'
            field= 'difference'

            diff_field = rfile.get_field(sweep, 'difference', copy=False)
            #  diff_field = Data['P_Radar'].rfile.fields['difference']['data']
            vminb, vmaxb = np.min(diff_field), np.max(diff_field)

        ###
        elif mom == 'Meso':
            c_label, c_scale = 'Azi Shear 'r' (s$^{-1}$)'+' '+Rname+' '+str(Scan_time), 'RdBu'
            field, vminb, vmaxb= 'Meso_azi', -10, 10 
        elif mom == 'Meso_test':
            c_label, c_scale = 'Azi Shear 'r' (s$^{-1}$)'+' '+Rname+' '+str(Scan_time), 'RdBu'
            field, vminb, vmaxb= 'Meso_test', -10, 10 

        elif mom == 'Ray_angle':
            vminb, vmaxb=  -90., 90
            field= mom 
            c_scale='twilight_shifted'
            c_label= 'Angle of Ray rel to Surge (deg)' 

        elif mom == 'Surge_subset_'+str(Surge_ID):
            vminb, vmaxb= -30., 30.
            field= mom 
            c_scale = cmocean.cm.balance
            c_label = 'Vort smooth 'r' [$\frac{1}{s}$]'
        ###
        else: print("Hey what just happened!\n Check the Radar moments for spelling")
        ##############################################################################
        #  print("This is the id for "+ field+' :'+str(id(Data['P_Radar'].rfile.get_field(Data['P_Radar'].swp, field, copy=False))))
        #  print(Data['P_Radar'].rfile.info(level='compact'))
        ## set colors to highlight where the max/min are exceded
        colormap = copy.copy(matplotlib.cm.get_cmap(name=c_scale))
        if mom == 'Ray_angle':
            colormap.set_over('green')
        else:
            colormap.set_under('aqua')
            colormap.set_over('magenta')


        ## Plot the radar
        #################
        #  ax_n.set_title(p_title, y=-.067, fontdict=self.Radar_title_font)
        DISPLAY.plot_ppi_map(field, sweep, ax=ax_n, cmap=colormap, vmin=vminb, vmax=vmaxb,
                                  width=self.config.radar_controls['offsetkm'][self.config.Radar_Plot_Type]*2000,
                                  height=self.config.radar_controls['offsetkm'][self.config.Radar_Plot_Type]*2000,
                                  title_flag=False, colorbar_flag=False, embelish=False)
        if mom == 'Ray_angle':
            vminb, vmaxb= -30., 30.
            c_scale2 = cmocean.cm.balance
            colormap2 = copy.copy(matplotlib.cm.get_cmap(name=c_scale))
            colormap2.set_under('aqua')
            colormap2.set_over('magenta')
            DISPLAY.plot_ppi_map('Surge_subset_'+str(Surge_ID), sweep, ax=ax_n, cmap=colormap2, vmin=vminb, vmax=vmaxb,
                                      width=self.config.radar_controls['offsetkm'][self.config.Radar_Plot_Type]*2000,
                                      height=self.config.radar_controls['offsetkm'][self.config.Radar_Plot_Type]*2000,
                                      title_flag=False, colorbar_flag=False, embelish=False, alpha=.3)
        ## Plot Contours  
        if self.config.overlays['Contour']['Lines']==True:
            unsmoothed_contourdata = rfile.get_field(sweep, self.config.overlays['Contour']['Var'])
            # smooth out the lines
            contourdata = ndimage.gaussian_filter(unsmoothed_contourdata, sigma=2)

            lat, lon, _ = rfile.get_gate_lat_lon_alt(sweep)
            levels = np.arange(-40, 40, self.config.overlays['Contour']['Interval'])
            contours=ax_n.contour(lon, lat, contourdata, levels, linewidths=1, antialiased=True, alpha=1, colors='k',
                                  linestyles='solid', transform=self.Proj)
            if self.config.overlays['Contour']['Label']== True: 
                ax_n.clabel(contours, levels,  fmt='%r', inline=True, fontsize=self.Contour_label_font['fontsize'])
        
        ## PLOT SPATIAL INDICATORS 
        ##########################
        # Det the domain (spacial extent) to be plotted on this axes 
        if mom == 'zoom' or ZOOM == True: 
            ax_n.set_extent(self.Surges[Surge_ID]['zoom_Domain'])
        else:
            ax_n.set_extent(self.Domain)

        # if the zoom will be plotted elsewhere in the fig but the current axes is the 'full' domain show a box that outlines the zoom area
        if ('zoom' in self.config.r_mom[:]) and (mom != 'zoom') and (ZOOM == False):
            if self.config.overlays['surge_lines']['zoom_boxes']== True:
                for Surge_ID in self.surge_ids[:]:
                    zoom_Domain= self.Surges[Surge_ID]['zoom_Domain']
                    rect = mpatches.Rectangle((zoom_Domain.xmin, zoom_Domain.ymin), abs(zoom_Domain.xmax-zoom_Domain.xmin),
                                              abs(zoom_Domain.ymax-zoom_Domain.ymin), fill=False, linewidth=2, transform=self.Proj)
                    ax_n.add_patch(rect)
        
        # Plot range rings on zoomed image
        if self.config.radar_controls['distance_indicator']['range_rings']== True:
            if mom == 'zoom':
                RANGE=np.repeat([rfile.range['data']], np.shape(self.gate_x)[0], axis=0) 
                levels = np.arange(-20000, 20000, self.config.radar_controls['distance_indicator']['ring_int'])
                contours=ax_n.contour(self.gate_lon, self.gate_lat, RANGE, levels, linewidths=2, linestyle=':',antialiased=True, 
                                      alpha=.6, colors='k', linestyles=':', transform=self.Proj)
                ax_n.clabel(contours, levels, fmt='%r', inline=True, fontsize=self.Contour_label_font['fontsize'])
        
        # Plot square grid on non-zoomed image
        if self.config.radar_controls['distance_indicator']['square_grid']== True:
            if mom != 'zoom' and ZOOM == False:
                self.tick_grid_settings(ax=ax_n, radar=True, scan_time=None, interval=.5*1000) # interval=10*1000)
                ax_n.grid(True)

        ## PLOT SURGE FEATURES 
        #######################
        ## Plot points along the surge boundary  
        if self.config.overlays['surge_lines']['Marker']== True:
            if Surge_ID == None: pass
            else:
                if self.config.overlays['surge_lines']['zoom_only']== True and mom != 'zoom':
                    pass
                else:
                    ax_n.plot(self.Surges[Surge_ID]['Center_pnts']['x'],
                              self.Surges[Surge_ID]['Center_pnts']['y'],
                              marker='.',mec='k',mfc='k',markersize=10, zorder=9, color='yellow',
                              path_effects=[PE.withStroke(linewidth=6, foreground='k')])#, c='red')
                    for o_dist in self.config.Surge_controls['Surges']['offset_dist']:
                        for o_dir in ['Plus', 'Minus']:
                            if o_dir == 'Plus': Mark='^' 
                            elif o_dir == 'Minus': Mark='D' 
                            ax_n.plot(self.Surges[Surge_ID]['Offset_pnts'][o_dist][o_dir]['x'],
                                      self.Surges[Surge_ID]['Offset_pnts'][o_dist][o_dir]['y'],
                                      marker=Mark, mec='k', mfc='grey', ms=6, color='grey', linestyle='solid', lw=1, zorder=10, 
                                      path_effects=[PE.withStroke(linewidth=3, foreground='k')])#, c='red')
                    for pnt in range(len(self.Surges[Surge_ID]['point_labels'])):
                        plt.scatter(self.Surges[Surge_ID]['Center_pnts']['x'][pnt],
                                    self.Surges[Surge_ID]['Center_pnts']['y'][pnt],
                                    zorder=10, s=40, 
                                    path_effects=[PE.withStroke(linewidth=5, foreground='k')])#, c='red')

                ax_n.scatter(self.Surges[Surge_ID]['Surge_Subset']['Max']['x'], self.Surges[Surge_ID]['Surge_Subset']['Max']['y'])
                ax_n.scatter(self.Surges[Surge_ID]['Surge_Subset']['Min']['x'], self.Surges[Surge_ID]['Surge_Subset']['Min']['y'])
                """Plot a line from slope and intercept"""
                x_vals = np.array(ax_n.get_xlim())
                y_vals = self.Surges[Surge_ID]['line']['intercept'] + self.Surges[Surge_ID]['line']['slope'] * x_vals
                print(y_vals)
                plt.plot(x_vals, y_vals, '--', linewidth=2, color='k') #,transform=self.Proj)
        
        if self.config.overlays['Mesos']['Marker']== True:
            if Meso_ID == None: pass
            else:
                for Meso in Meso_ID:
                    ax_n.plot(self.Mesos[Meso]['Center_pnt']['x'],
                              self.Mesos[Meso]['Center_pnt']['y'],
                              marker='.',mec='k',mfc='k',markersize=2, zorder=9, color='yellow',
                              path_effects=[PE.withStroke(linewidth=6, foreground='k')])#, c='red')
                    ax_n.scatter(self.Mesos['Max']['x'], self.Mesos['Max']['y'])
                    ax_n.scatter(self.Mesos['Min']['x'], self.Mesos['Min']['y'])

                    plt.plot(*self.Mesos[Meso]['Center_pnt']['Ring'].exterior.xy)
            ## PLOT PLATFORMS AS OVERLAYS(ie marker,colorline etc) ON RADAR
        ###############################################################
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
            #########
            '''Plot the filled line that indicates the thermo values on the Rplt'''
            #grab the subset of data of +- interval around radar scan
            p_sub, p_deploy = p.grab_pform_subset(p, Data, time_offset=self.config.overlays['Colorline']['cline_extent'])
            
            #if there is data for the platform that falls within the time and location of interest
            if p_deploy == True:
                ##READ in data from p_sub
                x, y = p_sub['lon'].values, p_sub['lat'].values

                if self.config.overlays['Colorline']['Pert'] == True: 
                    C = p_sub[self.config.overlays['Colorline']['Var']+'_pert'].values 
                    CSCALE='RdBu_r'
                else: 
                    C = p_sub[self.config.overlays['Colorline']['Var']].values
                    CSCALE=cmocean.cm.thermal

                ##PLOT the line that extends +/- min from the platform location;
                if self.config.overlays['Colorline'][p.type] == True:
                    #  The colorfill indicates values of the specifiec p_var (ie Thetae etc)
                    ax_n.scatter(x, y, c=C, cmap=CSCALE, norm=TVARS.norm, alpha=.7, marker='o',
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
                C_Point=p_sub.loc[p_sub['datetime'].sub(Scan_time).abs().idxmin()]
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
                    lat2, lon2 = Platform.getLocation(pform, radius, scan_time=Scan_time, given_bearing=bearing)
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
        
        if self.config.overlays['Tracks']['Marker_test']== True:
            track_type =['pnts', 'polys', 'paths'] 
            track_data = self.config.g_TORUS_directory+day+'/data/reports/extractDamage/nws_dat_damage_'+ track_type

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
            #  print(shapeRecs[0].record) #  print(shapeRecs) #  print(len(shapeRecs))  #  print(shapeRecs.record)
            #  print('999999999')
            # quick example plot to make sure data pulled correctly 
            for tor in np.arange(0, len(tor_tracks)):
                ax_n.scatter(tor_tracks.loc[tor]['lons'], tor_tracks.loc[tor]['lats'], label=tor_tracks.loc[tor]['rating'], transform=self.Proj)
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
                        #  xs.append(temp[0]) #  ys.append(temp[1])
                    #  plt.plot(np.asarray(xs)/1000, np.asarray(ys)/1000, color=colors[icol], label=tor_event)
                #  icol+=spyi1
        #  print('HHHHHHHHHHHHHHHHH')

        # * * *
        ## PLOT BACKGROUND FEATURES
        if self.config.overlays['GEO']['OX'] == True: P_Radar.plot_topo(Master_Plt, ax_n)

        # * * *
        ## DEAL WITH COLORBARS
        # Attach colorbar to each subplot
        if self.config.radar_controls['colorbars']['mom'] == True:
            divider = make_axes_locatable(plt.gca())
            c_ax = divider.append_axes("bottom", "5%", pad="2%", axes_class=plt.Axes)
            #  sm = plt.cm.ScalarMappable(cmap=c_scale, norm=matplotlib.colors.Normalize(vmin=vminb, vmax=vmaxb))
            sm = plt.cm.ScalarMappable(cmap=colormap, norm=matplotlib.colors.Normalize(vmin=vminb, vmax=vmaxb))
            sm._A = []
            cb = plt.colorbar(sm, cax=c_ax, orientation='horizontal', extend='both', drawedges=testing_testing, boundaries=testing_bound)
            cb.set_label(label=c_label, size=self.Radar_title_font['fontsize']) #, weight='bold')

        if self.config.radar_controls['colorbars']['thermo'] == True and thermo_cbar == True:
            #set up the colorbar for the thermo colorline overlay
            #  divider = make_axes_locatable(plt.gca())
            c_ax2 = divider.append_axes("right", "5%", pad="2%", axes_class=plt.Axes)
            Thermo_cbar = plt.colorbar(TVARS.CS3, cax=c_ax2, label=TVARS.Tvar_lab, ticks=MaxNLocator(integer=True))#use_gridspec=True)

        # * * *
        ## SET UP LEGENDS
        if leg == True: #this means you are currently making the left subplot
            #add legend for platform markers
            l = ax_n.legend(loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(-0.07,.5),
                            markerscale=1.1, labelspacing=.7, handletextpad=1.1, handlelength=.1)#, title="Platforms")
            l.set_title("Platforms", prop=self.leg_title_font)
            #  l._in_layout=False
            #  print(vars(l))

        ax_n.label_outer()
        ###################
        if self.config.print_long == True: print('~~~~~~~~~~~Made it through radar_subplots~~~~~~~~~~~')

    # * * * * * * * * * * * * * * * * * * * * * * * *
    # * * * * * * * * * * * * * * * * * * * * * * * *
    def line_plot(self, ts, ax_n, fig, TVARS, Data, Surge_ID, testing=False, deriv=False):
        ''' Plot a time series (sub)plot; could be of pvar info from various intstruments or of wind info
        ----
        INPUTS: ts: ........fill in
                Data: dict as described in plotting defn
                print_long & e_test: bool str as described in the ppi defn '''
        if self.config.print_long == True: print('~~~~~~~~~~~Made it into line_plot~~~~~~~~~~~~~~~~~')
        #  t_TS.T_Plt_settings(ts,ax_n)
        TSleg_elements = [] #empty list to append the legend entries to for each subplot that is actually plotted
        ## MAKE THE TIMESERIES
        #### * * * * * * * * *
        if ts == 'histogram':
            hist, bins = np.histogram(rfile.fields[self.config.lineplt_control['H']['var']]['data'][~np.isnan(
                rfile.fields[self.config.lineplt_control['H']['var']]['data'])], bins=150)
            bins = (bins[1:]+bins[:-1])/2.0
            ax_n.plot(bins, hist)
            ax_n.set_xlim(self.config.lineplt_control['H']['xmin'], self.config.lineplt_control['H']['xmax'])
            #  ax_n.set_yscale('log')
            #  ax_n.axvline(x=3,c='k')
            ax_n.set_xlabel(self.config.lineplt_control['H']['var'])
            ax_n.set_ylabel('Count')

        elif ts == 'clicker':
            WL=19
            if testing ==False:
                #  rmom=Data['P_Radar'].rfile.fields['Surge_subset_'+str(Surge_ID)]['data']
                rmom=self.Surges[Surge_ID]['Surge_Subset']['R_data']
                range_bins_alongray = Data['P_Radar'].rfile.range['data']
                print(np.shape(range_bins_alongray))
                print(np.shape(rmom))

                for ray in self.Surges[Surge_ID]['Surge_Subset']['Valid_Rays']:
                    #  r_bin_sub.append(range_bins_alongray[ranges_ind[pnt]-Max_obins:ranges_ind[pnt]+Max_obins])
                    data_subset = scipy.signal.savgol_filter(rmom[ray,:], window_length=WL, polyorder=2)#, deriv=1)#, axis=A)
                    print(data_subset)
                    print(np.shape(data_subset))
                    plt.plot(range_bins_alongray, data_subset)#, label=self.Surges[Surge_ID]['point_labels'][pnt])
                    #  plt.plot(self.Surges[Surge_ID]['Subsets']['Range_bins'][ray,:], data_subset)#, label=self.Surges[Surge_ID]['point_labels'][pnt])
                '''
                for pnt in range(len(self.Surges[Surge_ID]['point_labels'])):
                    self.Surges[Surge_ID]['Subsets']['Radar_data'][pnt,:]
                    data_subset = scipy.signal.savgol_filter(self.Surges[Surge_ID]['Subsets']['Radar_data'][pnt,:], window_length=WL, polyorder=2)#, deriv=1)#, axis=A)
                    plt.plot(self.Surges[Surge_ID]['Subsets']['Range_bins'][pnt,:], data_subset)#, label=self.Surges[Surge_ID]['point_labels'][pnt])
                '''
                '''
                for pnt in range(len(self.Surges[Surge_ID]['point_labels'])):
                    data_subset = scipy.signal.savgol_filter(self.Surges[Surge_ID]['Subsets']['Radar_data'][pnt,:], window_length=WL, polyorder=2)#, deriv=1)#, axis=A)
                    plt.plot(self.Surges[Surge_ID]['Subsets']['Range_bins'][pnt,:], data_subset)#, label=self.Surges[Surge_ID]['point_labels'][pnt])
                    
                    plt.scatter(self.Surges[Surge_ID]['Center_pnts']['range_value'][pnt],
                                self.Surges[Surge_ID]['Center_pnts']['radar_value'][pnt],
                                s=40)#c='red')
                for o_dist in self.config.Surge_controls['Surges']['offset_dist']:
                    for o_dir in ['Plus', 'Minus']:
                        if o_dir == 'Plus': Mark='^' 
                        elif o_dir == 'Minus': Mark='D' 
                        plt.scatter(self.Surges[Surge_ID]['Offset_pnts'][o_dist][o_dir]['range_value'],
                                    self.Surges[Surge_ID]['Offset_pnts'][o_dist][o_dir]['radar_value'],
                                    s=40, marker=Mark, c='grey')
                '''
            else:
                for pnt in range(len(self.Surges[Surge_ID]['point_labels'])):
                    data_subset = scipy.signal.savgol_filter(self.Surges[Surge_ID]['Subsets']['Radar_data'][pnt,:], window_length=WL, polyorder=2, deriv=1)#, axis=A)
                    #  grad_data_subset=np.gradient(data_subset)
                    plt.plot(self.Surges[Surge_ID]['Subsets']['Range_bins'][pnt,:], data_subset)#, label=self.Surges[Surge_ID]['point_labels'][pnt])

            #  plt.legend()
            ax_n.grid(axis='both')
            ax_n.tick_params(which='both', axis='y', direction='in')
            ax_n.tick_params(which='major', axis='both', direction='in', width=2, length=7, color='grey')
            ax_n.tick_params(which='both', axis='y', labelrotation=90)
            #  ax.tick_params(which='minor', axis='both', direction='in', width=2, length=7, color='grey')
            plt.xlabel('Range from Radar (along ray) (m)')
            #  ax_n.yaxis._in_layout=False
            #  print(vars(ax_n.yaxis._in_layout=False)
            if testing ==False:
                plt.ylabel('Radial Velocity (m/s)')
            else:
                #  ax_n.set_ylim(-1,1)
                plt.ylabel('Radial Gradient Velocity (m/s)')
            ax_n.axhline(0, color='grey', linewidth=3, alpha=.5, zorder=1)
            ax_n.invert_xaxis()


                
        elif ts in ['Thetae', 'Thetav', 'Wind']:
            if ts in ['Thetav', 'Thetae']:
                for p in Data.values():
                    if isinstance(p, Torus_Insitu):
                        if self.config.print_long == True: print('Plotting '+str(p.name)+' on time series')
                        ## If the platform matches a type listed for masking it uses masked data for the time series; 
                        #  else use the unmasked data
                        if p.type in self.config.lineplt_control['Mask'][:]:
                            if ts == "Thetae":   plotting_data= p.Thetae_ts_mask_df
                            elif ts == "Thetav":   plotting_data= p.Thetav_ts_mask_df
                        else: plotting_data= p.df[ts].values
                        
                        if self.config.overlays['Colorline']['Pert'] == True: 
                            var=ts+'_pert'
                            plotting_data=p.df[var].values
                            #  base_state =get_timeave_previous(p, Data, 60)
                            #  plotting_data= plotting_data[:] - base_state
                            ax_n.axhline(0, color='grey', linewidth=3, alpha=.5, zorder=1)
                        
                        if deriv == True:  
                            if self.config.overlays['Colorline']['Pert'] == True: 
                                var=ts+'_pert_der'
                            else:
                                var=ts+'_der'
                            plotting_data=p.df[var].values
                            ax_n.axhline(0, color='grey', linewidth=3, alpha=.5, zorder=1)
                        
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
                p = Data[self.config.lineplt_control['Wind_Pform']]
                if self.config.print_long == True: print('Plotting '+str(p.name)+' on time series')

                ax_n.plot(p.df['datetime'], p.df['spd'], label= 'Wind Spd')
                ax_n.fill_between(p.df['datetime'], p.df['spd'], 0)
                ax_n.set_ylabel('Wind Speed (kn)')
                self.tick_grid_settings(ax=ax_n, ts=ts,scan_time=P_Radar.Scan_time)
                TSleg_entry = Line2D([], [], label='Wind Spd', linewidth=12, color='tab:blue')
                TSleg_elements.append(TSleg_entry)

                ax_2 = ax_n.twinx()
                ax_2.plot(p.df['datetime'], p.df['dir'], '.k', linewidth=.05, label='Wind Dir')
                ax_2.set_ylabel('Wind Dir ($^{\circ}$)')
                self.tick_grid_settings(ax=ax_2, ts=ts, twinax=True, scan_time=P_Radar.Scan_time)
                TSleg_entry = Line2D([], [], marker='.', color='black', label='Wind Dir', markersize=26)
                TSleg_elements.append(TSleg_entry)

            ## Plot legend
            leg = ax_n.legend(handles= TSleg_elements, scatterpoints=3, loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(-0.04,.5))
            if ts == 'Wind':
                leg.set_title(self.config.lineplt_control['Wind_Pform'], prop=self.leg_title_font)
                leg.remove()
                ax_2.add_artist(leg)

            #if plotting more than one timeseries only include the xaxis label & ticks to the bottom timeseries
            ax_n.label_outer()

            #if you are not making the bottom timeseries
            if ax_n.get_subplotspec().rowspan.stop != len(self.config.Line_Plot):
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

            if deriv == True:
                leg.remove()
                ax_n.set_ylim(-.2,.2)
            ##Label Formatter
            # # # # # # # # #
            ax_n.set_xlabel('Time (UTC)')
            ax_n.yaxis.set_label_coords(-.03, .5)
            #  make_axes_area_auto_adjustable(ax_n)

            ###################
        if self.config.print_long == True: print('~~~~~~~~~~~Made it through line_plot~~~~~~~~~~~~~~')
    
    # * * * * * * * * * * * * * * * * * * * * * * * *
    # * * * * * * * * * * * * * * * * * * * * * * * *
    def clicker_stats(self, ax_n, fig, TVARS, Data, Surge_ID):
        if Surge_ID == None: pass
        else:
            ax_n.text(-2.4, .98, 'Surge ID: '+str(Surge_ID), transform=ax_n.transAxes, ha="left", va="top", weight='bold', size=self.leg_title_font['size'] )
            current_surge_df= self.valid_surges[self.valid_surges['Surge_ID'] == Surge_ID]
            rname=current_surge_df.loc[0,'RName']
            ax_n.text(-2.4, .88, 'Radar origin: '+str(rname)+', '+ str(current_surge_df.loc[0,'Tilt']), transform=ax_n.transAxes, ha="left", va="top")

            ax_n.text(-2.4, .78, 'Points Selected: '+str(len(self.Surges[Surge_ID]['point_labels'])), transform=ax_n.transAxes, ha="left", va="top")
            
            init_point = Point(self.Surges[Surge_ID]['Center_pnts']['x'][0], self.Surges[Surge_ID]['Center_pnts']['y'][0])
            final_point = Point(self.Surges[Surge_ID]['Center_pnts']['x'][-1], self.Surges[Surge_ID]['Center_pnts']['y'][-1])
            straight_dist=init_point.distance(final_point)
            ax_n.text(-2.4, .68, 'Straight length: '+str(round(straight_dist)), transform=ax_n.transAxes, ha="left", va="top")
            
            ax_n.text(-2.4, .58, 'length: '+str(round(self.Surges[Surge_ID]['Center_pnts']['line_object'].length)), transform=ax_n.transAxes, ha="left", va="top")
                #  line.boundary #  line.centroid #  line.centroid.wkt #  line.interpolate(.75, normalized=True).wkt
            
            ax_n.text(-2.4, .48, 'area: '+str(round(self.Surges[Surge_ID]['Polygon'].area)), transform=ax_n.transAxes, ha="left", va="top")
                #  poly.minimum_rotated_rectangle

            ax_n.text(-2.4, .38, 'Slope: '+str(round(self.Surges[Surge_ID]['line']['slope'],3)), transform=ax_n.transAxes, ha="left", va="top")
            ax_n.text(-2.4, .28, 'Intercept: '+str(round(self.Surges[Surge_ID]['line']['intercept'])), transform=ax_n.transAxes, ha="left", va="top")
            ax_n.text(-2.4, .18, 'Max Vel: '+str(round(self.Surges[Surge_ID]['Surge_Subset']['Max']['value'])), transform=ax_n.transAxes, ha="left", va="top")
            ax_n.text(-2.4, .08, 'Min Vel: '+str(round(self.Surges[Surge_ID]['Surge_Subset']['Min']['value'])), transform=ax_n.transAxes, ha="left", va="top")
            if self.num_of_mesos != 0:
                ax_n.text(-2.4, -.02, 'Dist to Meso Ring: '+str(round(self.Surges[Surge_ID]['Center_pnts']['line_object'].distance(self.Mesos[self.meso_ids[0]]['Center_pnt']['Ring']))),
                          transform=ax_n.transAxes, ha="left", va="top")
                ax_n.text(-2.4, -.12, 'Dist to Meso Point: '+str(round(self.Surges[Surge_ID]['Center_pnts']['line_object'].distance(self.Mesos[self.meso_ids[0]]['Center_pnt']['Point_object']))),
                         transform=ax_n.transAxes, ha="left", va="top")
                ax_n.text(-2.4, -.22, 'Max in Meso Ring: '+str(round(self.Mesos['Min']['value']))
                          ,transform=ax_n.transAxes, ha="left", va="top")
                ax_n.text(-2.4, -.32, 'Min in Meso Point: '+str(round(self.Mesos['Max']['value']))
                          ,transform=ax_n.transAxes, ha="left", va="top")
            ### 

            bb = Bbox([[-2.5, -.4], [.9, .99]])
            fancy = mpatches.FancyBboxPatch((bb.xmin, bb.ymin), bb.width, bb.height, boxstyle=mpatches.BoxStyle('Round', pad=.02), fc='white', clip_on=False, path_effects=[PE.withSimplePatchShadow(alpha=.8)])
            ax_n.add_patch(fancy)

            #  ax_n.set_xlabels(None)
            ax_n.tick_params(which='both', axis='both', left=False, bottom=False, labelleft=False, labelbottom=False)
            ax_n._frameon=False
                
            #  print(vars(ax_n))
                #  plt.plot(*line.xy)
                #  buf=line.buffer(5000, cap_style=3)
                #  plt.plot(*buf.exterior.xy)

    ##########################################################################
###########################
## PLOTTING DEFINTIONS  ###
###########################
def plotting(config, day, Data, TVARS, start_comptime, tilt=None):
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
    PLT = Master_Plt(day, Data, config, tilt)

    ## Make the Radar Subplots (if applicable)
    #### * * * * * * * * * * * * * * * * * ** * * * *
    if len(config.r_mom) != 0:
        if 'sim_vel' in config.r_mom:
            sd_test = singledop.SingleDoppler2D( L=30.0, radar=Data['P_Radar'].rfile, range_limits=[1, 6], filter_data=True,
                    filter_distance=5, thin_factor=[4, 12], sweep_number=Data['P_Radar'].swp, name_vr='velocity')

            display = singledop.AnalysisDisplay(sd_test)
            #  display = singledop.AnalysisDisplay(Data['P_Radar'].rfile.fields[sd_test])
            display.four_panel_plot(scale=400, legend=20, return_flag=False, thin=6,
                        levels=-30.0+2.0*np.arange(31), name_vr='vel_fix', name_dz='refl_fix', dz_cmap='pyart_HomeyerRainbow', xlim=[-21,21], ylim=[-21,21])
        else:
            mom_index=0
            surge_id_index=0
            for row in range(PLT.R_rows):
                for col in range(PLT.R_cols):
                    try:
                        mom=config.r_mom[mom_index]
                        if PLT.num_of_surges== 0: Surge_ID = None 
                        else: Surge_ID = PLT.surge_ids[surge_id_index]
                    except:
                        if mom != 'Meso':
                            mom = 'pass'
                    
                    if mom != 'Meso':
                        RADAR= Data['P_Radar']
                    else: 
                        if PLT.num_of_mesos != 0: 
                            #  RNAME=PLT.valid_mesos['RName'][0]
                            #  print(PLT.valid_mesos)
                            #  RNAME=PLT.valid_mesos.RName.unique()[0]
                            #  TIME=PLT.valid_mesos.loc[0,'Scan_time'][0]
                            #  TIME=PLT.valid_mesos.Scan_time.unique()[0]
                            #  TIME=TIME.loc[0]
                            #  print(TIME)
                            #  TIME=dt.datetime.strptime(TIME, '%Y-%m-%d %H:%M:%S.%f')
                            #  WSR_info = det_nearest_WSR(Data[config.radar_controls['Centered_Pform']].df)
                             
                            #  radar_file = get_WSR_from_AWS(config, day, (TIME-timedelta(minutes=1)), (TIME+timedelta(minutes=1)), 'KLBB')
                            #  radar = read_from_radar_file(radar_file[0], True)
                            #  radar = radar_fields_prep(config, radar, 'WSR', 1, moment='Meso_azi')

                            #  print(radar.info(level='compact'))
                            #  RADAR={'name':RNAME, 'Scan_time':TIME, 'rfile':radar, 'swp':[0,1]}
                            RADAR=PLT.Mesos['RADAR']
                        else:
                            mom = 'pass'

                    print('Radar plot: '+ mom +', Outer GridSpec Pos: [0, :], SubGridSpec Pos: ['+str(row)+', '+str(col)+']')
                    #  print('MESO')
                    #  print(PLT.valid_mesos)
                    #  print('SURGES')
                    #  print(PLT.valid_surges)
                    if mom != 'pass':
                        ax_n= PLT.fig.add_subplot(PLT.r_gs[row, col], projection= PLT.R_Proj)
                        if row==0 and col==0: leg = True
                        else: leg = False
                        
                        if col == PLT.R_cols-1:  thermo_cbar = True
                        else: thermo_cbar = False

                        PLT.radar_subplots(RADAR, mom, day, ax_n, PLT.fig, TVARS, Data, leg, thermo_cbar, Surge_ID, PLT.meso_ids)

                    if mom == 'zoom': 
                        surge_id_index=surge_id_index+1
                        if surge_id_index >= (len(PLT.surge_ids)):  mom_index =mom_index+1
                    else:
                        mom_index =mom_index+1

                        

    ## Make the Times Series Subplots (if applicable)
    #### * * * * * * * * * * * * * * * * * ** * * * *
    if len(config.Line_Plot) != 0:
        #  for (subrow,), ts in np.ndenumerate(config.Line_Plot):
        lineplot_index=0
        surge_id_index=0
        for row in range(PLT.L_Rows):
            #  for col in range(PLT.L_Cols):
            print('Line plot: '+ config.Line_Plot[lineplot_index]+', Outer GridSpec Pos: [1, :], SubGridSpec Pos: ['+str(row)+', :]')
            if PLT.num_of_surges== 0: Surge_ID = None 
            else: Surge_ID = PLT.surge_ids[surge_id_index]

            if row != 0 and config.Line_Plot[lineplot_index] != 'clicker': 
                ax_n=PLT.fig.add_subplot(PLT.l_gs[row,1], sharx=ax_n)
            else:
                if 'clicker' in config.Line_Plot and Surge_ID ==None: 
                    Make_Line_plt= False
                else: 
                    Make_Line_plt= True 
                    if config.Line_Plot[lineplot_index] == 'clicker': col =1
                    else: col= 0 
                    

            if Make_Line_plt == True: 
                ax_n= PLT.fig.add_subplot(PLT.l_gs[row, col])
                PLT.line_plot(config.Line_Plot[lineplot_index], ax_n, PLT.fig, TVARS, Data, Surge_ID)
                if config.lineplt_control['Deriv'] == True:  
                    print(PLT.L_Cols)
                    ax_n= PLT.fig.add_subplot(PLT.l_gs[row, 1])
                    PLT.line_plot(config.Line_Plot[lineplot_index], ax_n, PLT.fig, TVARS, Data, Surge_ID, deriv=True)

                if 'clicker' in config.Line_Plot: 
                    ax_n=PLT.fig.add_subplot(PLT.l_gs[row,2])
                    PLT.line_plot(config.Line_Plot[lineplot_index], ax_n, PLT.fig, TVARS, Data, Surge_ID, testing=True)

                    ax_n=PLT.fig.add_subplot(PLT.l_gs[row,0])
                    PLT.clicker_stats(ax_n, PLT.fig, TVARS, Data, Surge_ID)

                    ax_n= PLT.fig.add_subplot(PLT.l_gs[row, 3], projection= PLT.R_Proj)
                    leg=False
                    PLT.radar_subplots(Data['P_Radar'], 'Ray_angle', day, ax_n, PLT.fig, TVARS, Data, leg, False, Surge_ID, PLT.meso_ids, ZOOM=True)
                    #  PLT.radar_subplots('Surge_subset_'+str(Surge_ID), day, ax_n, PLT.fig, TVARS, Data, leg, Surge_ID, ZOOM=True)

                    surge_id_index=surge_id_index+1
                    if surge_id_index >=(len(PLT.surge_ids)): 
                        lineplot_index= lineplot_index + 1
                else: 
                    lineplot_index= lineplot_index + 1

    ## Finish plot
    #### * * * * *
    ### Get outdir and outname (aka file path and file name) for the correct image setup
    #if radar is included in the image
    if len(config.r_mom) != 0: 
        plt.suptitle(Data['P_Radar'].site_name+' '+str(tilt)+r'$^{\circ}$ PPI '+Data['P_Radar'].fancy_date_str, y=PLT.title_spacer, ha='center')
        #if both radar and timeseries are included in the image
        if len(config.Line_Plot) != 0:
            file_string = '_'.join(config.Line_Plot)
            if 'clicker' in config.Line_Plot:
                plt_type = 'clicker'
            elif 'histogram' in config.Line_Plot:
                plt_type = 'hist'
            else: 
                if ('Thetae' in config.Line_Plot) and ('Thetav' in config.Line_Plot): plt_type = 'both_Thetas'
                elif (len(config.Line_Plot) != 1) and ('Wind' in config.Line_Plot): plt_type = config.overlays['Colorline']['Var']+'/Wind'
                else: plt_type = config.overlays['Colorline']['Var']

                if config.overlays['Colorline']['Pert'] == True:
                    plt_type = plt_type+str('_pert')

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
    elif len(config.Line_Plot) != 0 and len(config.r_mom) == 0:
        plt.suptitle(day+' Time Series', y=title_spacer)
        
        #add the file name
        file_string = '_'.join(config.Line_Plot)
        if config.ts_extent ==None: file_ender='_full.png'
        else: file_ender='.png'
        output_name= day+'_TS_'+file_string+file_ender
        #path to file 
        plt_dir='TSeries/'
    
    
    if config.Surge_controls['Feature_IDing']['Make_new_pnts']['activate_clicker'] == True:
        if config.Surge_controls['Feature_IDing']['Make_new_pnts']['Type'] == 'Surge':
            csv_name= config.g_TORUS_directory+day+'/data/'+day+'_surge_pnts.csv'
            Does_csv_exist=os.path.isfile(csv_name)
            surge_name= 'NA'
            #  plt.show()
            select_points=input('Are you interested in selecting points for this time (Y/N):\n')
            if select_points == 'Y':
                AXES=PLT.fig.get_axes()
                while surge_name != 'DONE':
                    print('####')
                    surge_name=input('Type name of surge you are adding points for (1,2,etc):\n If you are done adding points type DONE:\n')
                    print('####')
                    
                    if surge_name== 'DONE':
                        break
                    else:
                        pts=clicker_defn(PLT, AXES)
                        make_csv(config, Data, pts, surge_name, tilt, csv_name, 'Surge', Does_csv_exist)
        elif config.Surge_controls['Feature_IDing']['Make_new_pnts']['Type'] == 'Meso':
            csv_name= config.g_TORUS_directory+day+'/data/'+day+'_meso_pnts.csv'
            Does_csv_exist=os.path.isfile(csv_name)
            meso_name= 'NA'
            select_points=input('Are you interested in selecting points for this time (Y/N):\n')
            if select_points == 'Y':
                plt.setp(plt.gca(), autoscale_on=False)
                AXES=PLT.fig.get_axes()
                while meso_name != 'DONE':
                    print('####')
                    meso_name=input('Type name of meso you are identifing (M1,M2,etc):\n If you are done adding points type DONE:\n')
                    print('####')
                    
                    if meso_name== 'DONE':
                        break
                    else:
                        pts=clicker_defn(PLT, AXES)
                        make_csv(config, Data, pts, meso_name, tilt, csv_name, 'Meso', Does_csv_exist)
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
            if config.radar_controls['Use_downloaded_files']== True:
                path = config.g_TORUS_directory + day+'/data/radar/Nexrad/Nexrad_files/dealiased_K*'
                r_files_path = sorted(glob.glob(path))
            else:
                r_files_path = det_nearest_WSR(Data[config.radar_controls['Centered_Pform']].df)
        return r_files_path
    
    ####
    def plot_radar_file(config, r_file, day, Data, TVARS, subset_pnames):
        print("\n*******************************\n plot_radar_file: scan file_name = {}\n".format(r_file))
        start_comptime = time.time()
        if config.Radar_Plot_Type == 'WSR_Plotting':
            is_WSR= True
            if config.radar_controls['Use_downloaded_files']== True:
                time_string= str(os.path.split(r_file)[1][14:29])
                rtime = datetime.strptime(time_string, "%Y%m%d_%H%M%S")
            else:
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
                radar = read_from_radar_file(config, r_file, is_WSR)
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
                plotting(config, day, Data, TVARS, start_comptime, tilt)
                print("done in plot_radar_file")


    #########
    if config.Radar_Plot_Type == 'WSR_Plotting':
        if config.radar_controls['Use_downloaded_files']== True:
            radar_files_paths = find_radar_files(config, day, Data)
        else:
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
    elif len(plot_config.Line_Plot) != 0:
        # Create Timeseries only plot ##
        #Only plot timeseries (this code isn't fully fleshed out but in theroy this code is built in such a way to allow for this)
        plot_time_series(plot_config, day, Data, TVARS)

###########################
plt.rcdefaults()
timer(totalcompT_start, time.time(), total_runtime=True)
print("ALL FINISHED")

tracemalloc.stop()

