#!/usr/bin/env python
# -*- coding: utf-8 -*-

#hellow
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
import fnmatch

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
import config #this is the file with the plotting controls to access any of the vars in that file use config.var
if config.country_roads == True: import osmnx as ox

print_long =config.print_long
e_test = config.e_test

## Read in defns that I have stored in another file (for ease of use/consistancy accross multiple scripts)
from read_pforms import pform_names, Add_to_DATA, Platform, Torus_Insitu, Radar, Stationary_Insitu
from read_pforms import error_printing, timer, time_in_range

totalcompT_start = time.time()
###########
# Classes
##########
class Thermo_Plt_Vars:
    def __init__ (self, Data):
        self.Te_lab, self.Tv_lab = "Equi. Pot. Temp [K]", "Vir. Pot. Temp [K]"
        self.Tv_GMin, self.Tv_GMax = self.find_global_max_min('Thetav', Data)
        self.Te_GMin, self.Te_GMax = self.find_global_max_min('Thetae', Data)

        #settings for the thermo var being plotted on radar subplots
        if len(config.r_mom) != 0:
            if config.R_Tvar == "Thetae":
                lab, GMin, GMax = self.Te_lab, self.Te_GMin, self.Te_GMax
            elif config.R_Tvar == "Thetav":
                lab, GMin, GMax = self.Tv_lab, self.Tv_GMin, self.Tv_GMax

            self.R_Tvar_lab, self.R_Tvar_GMin, self.R_Tvar_GMax = lab, GMin, GMax

            ## Make a dummy plot to allow for ploting of colorbar
            #  cmap=plt.get_cmap('rainbow')
            self.pvar_color_setup(colorbar=True)

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
        if len(val_hold) != 0: global_min, global_max = min(val_hold), max(val_hold)
        return global_min, global_max

    # * * *
    def pvar_color_setup(self, cmap=cmocean.cm.thermal, df=None, colorbar=False, colorramp=False):
        if colorbar==True:
            Z = [[0,0],[0,0]]
            levels = np.arange(self.R_Tvar_GMin, self.R_Tvar_GMax+1, 1)

            self.CS3 = plt.contourf(Z, levels, cmap=cmap)
            plt.clf()

        if colorramp==True:
            self.C = cmap((df[config.R_Tvar].values - self.R_Tvar_GMin) / (self.R_Tvar_GMax - self.R_Tvar_GMin))
            return self.C

    # * * *
    def det_TS_Tvar(self, var):
        if var == "Thetae":
            lab, GMin, GMax = self.Te_lab, self.Te_GMin, self.Te_GMax
        elif var == "Thetav":
            lab, GMin, GMax = self.Tv_lab, self.Tv_GMin, self.Tv_GMax
        self.TS_Tvar_lab, self.TS_Tvar_GMin, self.TS_Tvar_GMax = lab, GMin, GMax

#############################################################################################
#####################################################
# Read in Data that does not change for each image ##
#####################################################
print('\nRead in Torus Platforms')
#print(pform_names('ALL')) #list of ALL possible platforms
subset_pnames = [] #array of platforms that actually have data for the time range we are plotting
Data = {} #dictionary in which all the class objects will be stored (will contain platform data, locations, plotting variables etc)

## Read in the data for the TORUS Insitu platforms (if available)
Data, subset_pnames = Add_to_DATA('TInsitu', Data, subset_pnames, config.print_long)
## Establish info for the plotting variable (aka max min etc)
TVARS=Thermo_Plt_Vars(Data)

#  print('\nRead in Stationary Platforms Arrays')
## Read in the data for the Stationary Array platforms (if available)
#  Data, subset_pnames = Add_to_DATA('STN_I', Data, subset_pnames, config.print_long)



def grab_sub(pform, lower_tbound, upper_tbound):
    #  lower_tbound, upper_tbound = (pform.Scan_time-dt.timedelta(minutes=time_offset)), (pform.Scan_time+dt.timedelta(minutes=time_offset))
    df_sub = pform.df.loc[(pform.df['datetime'] >= lower_tbound) & (pform.df['datetime'] <= upper_tbound)]
    return df_sub

 

day= '20190608'
tstart = dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), 20, 0, 0)
tend = dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),20,30,0)
truth_probe = 'FFld' 
print(tstart)
print(type(tstart))

            
thetae_df = pd.DataFrame()
thetav_df = pd.DataFrame()
'''
        for pname in pform_names('TInsitu'):
            #
            if pname in pform_names('NSSL'):
                if config.NSSLm == True: read_in_data= True
                else: read_in_data= False
            if pname in pform_names('UNL'):
                if config.NEBm == True: read_in_data= True
                else: read_in_data= False
            if pname == 'UAS':
                if config.UASm == True: read_in_data= True
                else: read_in_data= False

            # if you do want to read in data
            if read_in_data == True:
                #for each platform test to see if we have data for that day
                data_avail = Platform.test_data(pname)
                if data_avail == True:
                    subset_pnames.append(pname) #append the pname to the subset_pnames list
                    #  load data for the pform (aka initialize an object of the appropriate class); place in dict with key of pname
                    Data.update({pname: Torus_Insitu(pname)})
                    if print_long == True: print("Data can be read in for platform %s" %(pname))
                else:
                    if print_long == True: print("No data available to be read in for platform %s" %(pname))
            elif read_in_data == False:
                    if print_long == True: print("Did not attempt to read in data for platform %s" %(pname))

'''
Data_sub={}
print(Data)
#comparision probe
#grab the subset of data of +- interval around radar scan
#  truth_sub = grab_sub(Data[truth_probe], tend, tstart)
#  print(truth_sub)
for pname in pform_names('TInsitu'):
    for p in Data.values():
        if pname == p.name:
            print(pname)
            print(Data[pname].df)
            comp_sub = grab_sub(Data[pname], tend, tstart)
            print(comp_sub)
            Data_sub.update({pname:comp_sub}) 

print(Data_sub)
'''
for p in Data.values():
    if isinstance(p, Torus_Insitu):
        print(p.name)
        print(vars(p))
        comp_sub = grab_sub(p, tend, tstart)
        print(comp_sub.Thetav)
        thetae_df = thetae_df.assign(name = comp_sub.Thetae)
        thetav_df = thetav_df.assign(name = comp_sub.Thetav)
print('999999')
print(thetae_df)
print(thetav_df)
'''        
