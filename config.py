#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import needed modules
######################
import datetime as dt

import os
from pathlib import Path
from os.path import expanduser

# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## VARIABLES ##
###############
## File paths
# ************
#  Radar_Plot_Type = 'NOXP_Plotting'
#  Radar_Plot_Type = 'KA_Plotting'
Radar_Plot_Type = 'WSR_Plotting'
print_radar_info= False

nCPU=1

if Radar_Plot_Type == 'KA_Plotting':
    day = '20190517' #'YYYYMMDD'
    R_Tvar = "Thetae" #which var to plot (current options; Thetae, Thetav)
    offsetkm = 21 #21 is the best for KA
    Centered_Pform = 'Prb1' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    Wind_Pform = 'Prb1'
    p_tilt = 1 #what radar elevation tilt (deg) do you want to plot

if Radar_Plot_Type == 'NOXP_Plotting':
    day = '20190615' #'YYYYMMDD'
    R_Tvar = "Thetav" #which var to plot (current options; Thetae, Thetav)
    offsetkm = 37 #21 is the best for KA
    Centered_Pform = 'P_Radar' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    Wind_Pform = 'Prb1'
    p_tilt = 0.5 #what radar elevation tilt (deg) do you want to plot

if Radar_Plot_Type == 'WSR_Plotting':
    day = '20190527' #'YYYYMMDD'
    R_Tvar = "Thetav" #which var to plot (current options; Thetae, Thetav)
    offsetkm = 60 #21 is the best for KA
    Centered_Pform = 'Prb1' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    Wind_Pform = 'Prb1'
    p_tilt = 0.5 #what radar elevation tilt (deg) do you want to plot


## Timeframe
# ***********
#  day = '20190520' #'YYYYMMDD'
#crop the start or end time to a time you specify (set to None to do the whole day)
#  tstart = dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), 18, 30, 0)
#  tend = dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),22,45,0)
tstart, tend = None, None
#  tend = None

#g_mesonet_directory ='/Users/severe2/Research/TORUS_Data'
#g_plots_directory ='/Users/severe2/Research/TORUS_Data'
#g_download_directory ='/Users/severe2/Research/TORUS_Data'
#filesys='/Users/severe2/Research/'
#temploc='/Volumes/Samsung_T5/Research/TORUS_Data/'

g_home = expanduser("~")

#g_root = g_home + '/radar_data'   # Chris
#g_root = g_home + '/Research'     # Ellie Macbook
g_root = '/Volumes/Samsung_T5/Research' # Ellie External SSD

g_mesonet_directory = g_root + '/TORUS_Data/'   # Our Data
g_download_directory = g_root + '/downloads/'    # Download Storage
#  g_plots_directory = g_root + '/plots/' # Plot Output
g_plots_directory = g_mesonet_directory # Plot Output
g_cache_directory = g_root + '/cache/' # Cache
g_roads_directory = g_root + '/roads/' # Roads data

#
# Remember to erase all plots for completely new versions as they are not regenerated if the file already exists
# Plots are stored in g_root + /plots
#

# Setup Function Cache for speedup
if not os.path.exists(g_cache_directory): Path(g_cache_directory).mkdir(parents=True)


## Plot layout controls
# ***********************
#  (some of these options will need to be expanded upon to work.... but hypothetically the code is set up to allow these options(with some modifications))
#would you like to plot radar images? (set to False to only plot timeseries)
#would you like to plot the timeseries? (set to False to only plot radar)
r_plotting, t_plotting = True, True

## Plot visual controls (aka what is plotted)
# *******************************************
#  p_var = "Thetae" #which var to plot (current options; Thetae, Thetav)
#  p_var = "Thetav" #which var to plot (current options; Thetae, Thetav)
#  offsetkm = 30 #21 is the best for KA
#  offsetkm = 40 #21 is the best for KA
#  Centered_Pform = 'P_Radar' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
#  Centered_Pform = 'Prb1' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
#  Wind_Pform = 'Prb1'
#  p_tilt = .5 #what radar elevation tilt (deg) do you want to plot
#  p_tilt = .5 #what radar elevation tilt (deg) do you want to plot
rhi_ring = True #do you want the rhi spokes for the KA radars
r_mom = ['refl', 'vel'] #list; the radar moments to plot on the left and right subplots respectively (current options are 'refl' or 'vel)

#which other platforms do you wish to plot as overlays to the radar image (aka as markers, colorlines etc) if data is available
KAm, WSRm, NOXPm, NSSLm, NEBm, UASm, ASOSm, MESONETSm = True, True, True, True, True, False, True, True
MESO_lab, WSR88D_lab, KA_lab, NOXP_lab, RHI_lab, TIn_lab, ASOS_lab = False, False, False, False, False, False, False
country_roads, hwys, county_lines, state_lines = False, False, False, False #background plot features
cline_extent = 10 #how long would you like the colorlines for each platforms to cover +-x min (also control width of greybox on tseries)
ts_extent = 30 #so actually 60 min
#  ts_extent = None #so actually 60 min

#list; include in the list the type names for the pforms you would like masked (this only controls whats plotted on the Tseries)
#  if none leave as blank list, will not work if a mask had not been previously applied to the dataset (which can be done in the classes)
Time_Series = ['Thetae', R_Tvar]
#  Time_Series = ['Wind', p_var]

TS_masking = ['NSSL']
NSSL_qcflags = ['qc1', 'qc2', 'qc3', 'qc4']
# would you like to restrict the time series to display only a timespan +/- x time from the scantime (functionalitly not coded yet)
#  pri_tseries='p_var' #what would you like to plot in the first time series (if num_of_tseries=1 this is the only tseries to be plotted)
#  sec_tseries='wind' #what would you like in the second time series (only used if num_of_tseries =2)


## Troubleshooting (True/False variable)
# *****************
#  Set to True to print statements when entering/leaving definitions
#  (helpful to understanding whats going on in the code but prints alot of lines to your terminal)
print_long = False
#  Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting
    ## there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = False
#  nCPU=-2
#  nCPU=1

