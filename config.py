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
## Plot layout controls
# ***********************
r_mom = ['refl', 'vel'] #list; the radar moments to plot on the left and right subplots respectively (current options are 'refl' or 'vel)
#  r_mom = []
#  Time_Series = ['Wind', R_Tvar]
Time_Series = ['Thetav']

## File paths
# ************
#  Radar_Plot_Type = 'NOXP_Plotting'
#  Radar_Plot_Type = 'KA_Plotting'
Radar_Plot_Type = 'WSR_Plotting'

nCPU = -2

if Radar_Plot_Type == 'KA_Plotting':
    day = '20190524' #'YYYYMMDD'
    R_Tvar = "Thetav" #which var to plot (current options; Thetae, Thetav)
    offsetkm = 21 #21 is the best for KA
    Centered_Pform = 'P_Radar' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    Wind_Pform = 'P_Radar'
    #  p_tilt = 1  #what radar elevation tilt (deg) do you want to plot
    p_tilt = 1.0  #what radar elevation tilt (deg) do you want to plot

if Radar_Plot_Type == 'NOXP_Plotting':
    day = '20190517' #'YYYYMMDD'
    R_Tvar = "Thetav" #which var to plot (current options; Thetae, Thetav)
    #  offsetkm = 37 #21 is the best for KA
    offsetkm = 60 #21 is the best for KA
    Centered_Pform = 'P_Radar' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    Wind_Pform = 'P_Radar'
    p_tilt = .5 #what radar elevation tilt (deg) do you want to plot

if Radar_Plot_Type == 'WSR_Plotting':
    day = '20190520' #'YYYYMMDD'
    R_Tvar = "Thetav" #which var to plot (current options; Thetae, Thetav)
    offsetkm =100#70 #21 is the best for KA
    Centered_Pform = 'CoMeT1' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    Wind_Pform = 'CoMeT1'
    p_tilt = 0.5 #what radar elevation tilt (deg) do you want to plot


## Timeframe
# ***********
#crop the start or end time to a time you specify (set to None to do the whole day)
#  tstart = dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), 21, 30, 0)
#  tend = dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), 23, 0, 0)
man_ylim= False


#g_root = g_home + '/radar_data'   # Chris
#g_root = g_home + '/Research'     # Ellie Macbook
g_home = expanduser("~")
g_root = '/Volumes/Samsung_T5/Research' # Ellie External SSD

g_TORUS_directory = g_root + '/TORUS_Data/'   # Our Data
g_csv_directory = g_root+ '/csv_files/'
g_cache_directory = g_root + '/cache/' # Cache
g_download_directory = g_cache_directory + 'downloads/' # Download Storage
g_roads_directory = g_root + '/roads/' # Roads data

# Setup Function Cache for speedup
if not os.path.exists(g_cache_directory): Path(g_cache_directory).mkdir(parents=True)

## Plot layout controls
# ***********************
#would you like to plot radar images? (set to False to only plot timeseries)
#  r_mom = ['refl', 'vel'] #list; the radar moments to plot on the left and right subplots respectively (current options are 'refl' or 'vel)
#  r_mom = []#list; the radar moments to plot on the left and right subplots respectively (current options are 'refl' or 'vel)

#list; include in the list the type names for the pforms you would like masked (this only controls whats plotted on the Tseries)
#  if none leave as blank list, will not work if a mask had not been previously applied to the dataset (which can be done in the classes)
#  Time_Series = ['Wind', R_Tvar]
#  Time_Series = [R_Tvar]

## Plot visual controls (aka what is plotted)
# *******************************************
#  p_var = "Thetae" #which var to plot (current options; Thetae, Thetav)
#  offsetkm = 30 #21 is the best for KA
#  Centered_Pform = 'P_Radar' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
#  Wind_Pform = 'Prb1'
#  p_tilt = .5 #what radar elevation tilt (deg) do you want to plot
rhi_ring = True #do you want the rhi spokes for the KA radars

#which other platforms do you wish to plot as overlays to the radar image (aka as markers, colorlines etc) if data is available
KAm, WSRm, NOXPm, NSSLm, NEBm, UASm, ASOSm, MESONETSm = True, True, True, True, True, False, True, True
MESO_lab, WSR88D_lab, KA_lab, NOXP_lab, RHI_lab, TIn_lab, ASOS_lab = False, False, False, False, False, False, False
country_roads, hwys, county_lines, state_lines = False, False, False, False #background plot features
cline_extent = 5 #how long would you like the colorlines for each platforms to cover +-x min (also control width of greybox on tseries)

# would you like to restrict the time series to display only a timespan +/- x time from the scantime 
ts_extent = 30 #so actually 60 min
#  ts_extent = None

TS_masking = ['NSSL']
NSSL_qcflags = ['qc1', 'qc2', 'qc3', 'qc4']

## Troubleshooting (True/False variable)
# *****************
#  Set to True to print statements when entering/leaving definitions
#  (helpful to understanding whats going on in the code but prints alot of lines to your terminal)
print_long = False
#  Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting
    ## there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = False
print_radar_info= False 
#  nCPU=-2 #1
