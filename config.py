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
#  r_mom = ['refl','vel'] #list; the radar moments to plot on the left and right subplots respectively (current options are 'refl' or 'vel)
#  r_mom = ['vel', 'vel_texture', 'vel_text_texture', 'vel_aray_texture' ] #list; the radar moments to plot on the left and right subplots respectively (current options are 'refl' or 'vel)
#  r_mom = ['refl','vel', 'sim_vel']
#  r_mom = ['refl','vel', 'vel_grad', 'vel_smooth']
r_mom = ['vel', 'vel_smooth', 'vel_grad', 'diff_dbz']
#  r_mom = ['refl', 'vel','vel_texture','vel_grad1_1']
#  r_mom = ['refl', 'vel','vel_texture','specw']
#  r_mom=[]
#  Time_Series = ['Wind', 'Thetav']
#  Time_Series = ['Wind', R_Tvar]
#  Time_Series = ['Thetav']
Time_Series = []

#  Radar_Plot_Type = 'NOXP_Plotting'
Radar_Plot_Type = 'KA_Plotting'
#  Radar_Plot_Type = 'WSR_Plotting'

nCPU = 1# -2

if Radar_Plot_Type == 'KA_Plotting':
    R_Tvar = "Thetae" #which var to plot (current options; Thetae, Thetav)
    offsetkm = 21 #21 is the best for KA
    Centered_Pform = 'P_Radar' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    Wind_Pform = 'P_Radar'

if Radar_Plot_Type == 'NOXP_Plotting':
    R_Tvar = "Thetav" #which var to plot (current options; Thetae, Thetav)
    #  offsetkm = 37 #21 is the best for KA
    offsetkm = 30 #21 is the best for KA
    Centered_Pform = 'Prb1' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    Wind_Pform = 'P_Radar'

if Radar_Plot_Type == 'WSR_Plotting':
    R_Tvar = "Thetav" #which var to plot (current options; Thetae, Thetav)
    offsetkm = 40 #21 is the best for KA
    Centered_Pform = 'Prb2' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    Wind_Pform = 'Prb2'

#  p_tilt = 1  #what radar elevation tilt (deg) do you want to plot
p_tilt=[1.0]#0.5]#, 1.0, 1.5]
## Timeframe
# ***********
#crop the start or end time to a time you specify (set to None to do the whole day)
#  day_list=['20190608']
day_list= [ '20190524']
#  day_list= [ '20190517']
#  day = day_list[0]
#  tstart = dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), 23, 50, 0)
#  tend = dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), 21, 30, 0)
#  day_list= ['20190517','20190518','20190520', '20190523', '20190524','20190525','20190526',
           #  '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190520', '20190523', '20190524','20190525','20190526',
           #  '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= [ '20190608', '20190611', '20190613', '20190615']
#  day_list= [ '20190524', '20190613', '20190615']


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

## Plot layout controls
# ***********************
#would you like to plot radar images? (set to False to only plot timeseries)
#  r_mom = ['refl', 'vel'] #list; the radar moments to plot on the left and right subplots respectively (current options are 'refl' or 'vel)

#list; include in the list the type names for the pforms you would like masked (this only controls whats plotted on the Tseries)
#  if none leave as blank list, will not work if a mask had not been previously applied to the dataset (which can be done in the classes)
#  Time_Series = ['Wind', R_Tvar]

## Plot visual controls (aka what is plotted)
# *******************************************
Read_Data = {'Meso': True, 'Radar':True}

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

# would you like to restrict the time series to display only a timespan +/- x time from the scantime 
#  ts_extent = 30 #so actually 60 min
ts_extent = None
cline_extent = 5 #how long would you like the colorlines for each platforms to cover +-x min (also control width of greybox on tseries)
man_ylim= False

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
#  print_radar_info= True
print_radar_info=False 
#  nCPU=-2 #1


