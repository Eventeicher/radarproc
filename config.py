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
#list; the radar moments to plot on the left and right subplots respectively (leave list empty to plot only timeseries) 
#  r_mom=[]
r_mom = ['refl','vel']
#  r_mom = ['refl','vel', 'refl_despeck', 'vel_despeck'] 
#  r_mom = ['vel', 'vel_texture', 'vel_text_texture', 'vel_aray_texture' ] 
#  r_mom = ['refl','vel', 'sim_vel']
#  r_mom = ['refl','vel', 'vel_grad', 'vel_smooth']
#  r_mom = ['vel', 'vel_smooth', 'vel_savgol_grad', 'vel_grad']
#  r_mom = ['vel', 'vel_unfixed', 'diff_dbz']
#  r_mom = ['refl', 'vel','vel_texture','vel_grad1_1']
#  r_mom = ['refl', 'vel','vel_texture','specw']

#list; the timeseries to plot (leave list empty to plot only radar) 
Time_Series = []
#  Time_Series = ['Wind', 'Thetav']
#  Time_Series = ['Wind', R_Tvar]
#  Time_Series = ['Thetav']

# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
#number of CPUs to use (recomend either 1 or -2)
nCPU = 1 #-2

# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
# what radar elevation tilt (deg) do you want to plot
p_tilt=[0.5, 1.0, 1.5]

## What type of radar do you want to plot 
#  Radar_Plot_Type = 'NOXP_Plotting'
Radar_Plot_Type = 'KA_Plotting'
#  Radar_Plot_Type = 'WSR_Plotting'

if Radar_Plot_Type == 'KA_Plotting':
    R_Tvar = "Thetae" #which var to plot on the colorline overlay (current options; Thetae, Thetav)
    offsetkm = 21 #size of the plot domain (times by 2): 21 is the best for KA
    Centered_Pform = 'P_Radar' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')

if Radar_Plot_Type == 'NOXP_Plotting':
    R_Tvar = "Thetav" 
    offsetkm = 30 
    Centered_Pform = 'Prb1' 

if Radar_Plot_Type == 'WSR_Plotting':
    R_Tvar = "Thetav" #which var to plot (current options; Thetae, Thetav)
    offsetkm = 40 #21 is the best for KA
    Centered_Pform = 'Prb2' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')


# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## Timeframe
# ***********
## list of days to plot
day_list= ['20190528']
#  day_list= ['20190517','20190518','20190520', '20190523', '20190524','20190525','20190526',
             #  '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
## Crop the start or end time to a time you specify (comment out to use datafrom the whole day)
#  tstart = dt.datetime(int(day_list[0][0:4]), int(day_list[0][4:6]), int(day_list[0][6:8])+1, 0, 0, 0)
#  tend = dt.datetime(int(day_list[0][0:4]), int(day_list[0][4:6]), int(day_list[0][6:8]), 21, 30, 0)


# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## File paths
# ************
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
if not os.path.exists(g_cache_directory): Path(g_cache_directory).mkdir(parents=True, exist_ok=True)

# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## Data controls  
# *******************************************
#whether to bother trying to read in this datatype at all (timesaver) 
Read_Data = {'Meso': True, 'Radar':True}

#list; include in the list the type names for the pforms you would like masked (this only controls whats plotted on the Tseries)
TS_masking = ['NSSL']
NSSL_qcflags = ['qc1', 'qc2', 'qc3', 'qc4']

# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## Plot visual controls (aka what is plotted)
# *******************************************
#what to plot on top of the radar image 
overlays = {'KA': {'Marker': False, 'Label': False, 'RHI_ring': False, 'RHI_Lab': False},
            'WSR': {'Marker': False, 'Label': False},
            'NOXP': {'Marker': False, 'Label': False},
            'MAINR':{'Marker':False, 'Label': False}, # always leave this one False
            'NSSL': {'Marker': False, 'Label': False, 'Colorline': False},
            'UNL': {'Marker': False, 'Label': False, 'Colorline': False},
            'UAS': {'Marker': False, 'Label': False},
            'LSR': {'Marker': False, 'Label': False},
            'Tracks': {'Marker': False, 'Label': False},
            'STN_I': {'Marker': False, 'Label': False},
            'GEO': { 'OX': False, 'country_roads': False, 'hwys': False, 'county_lines': False, 'state_lines': False}
           }

## Time Series controls
# would you like to restrict the time series to display only a timespan +/- x time from the scantime 
#  ts_extent = 30 #so actually 60 min 
ts_extent = None
#if you are plotting a wind time series what platform would you like it to be 
Wind_Pform = 'Prb1' 
#how long would you like the colorlines for each platforms to cover +-x min (also control width of greybox on tseries)
cline_extent = 5
#hardcode a y axis extent for the timeseries
man_ylim= False


# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## Troubleshooting (True/False variable)
# *****************
#  Set to True to print statements when entering/leaving definitions
#  (helpful to understanding whats going on in the code but prints alot of lines to your terminal)
print_long = False
#  Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting
    ## there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = False
#  print_radar_info= True
print_radar_info=False 


