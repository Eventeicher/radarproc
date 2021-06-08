#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import needed modules
######################
import datetime as dt
import os
from pathlib import Path
from os.path import expanduser

# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
#number of CPUs to use (recomend either 1 or -2)
#  nCPU = -2
nCPU = 1#-2

# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## VARIABLES ##
###############
## Plot layout controls
# ***********************
#list; the radar moments to plot on the left and right subplots respectively (leave list empty to plot only timeseries) 
    #  sim_vel, refl_despeck, vel_despeck, vel_unmasked, vel_unfixed, vel_test,
    #  vel_savgol_grad, vel_savgol, vel_texture, vel_texture_dea, specw, vel_savgol_der,
    #  vel_savgol_der2, vel_smooth, az_shear, vel_grad, vel_unfixed, diff_dbz, vel_grad1_1
#  r_mom=[]
#  r_mom = ['refl','vel', 'vel_texture']
#  r_mom = ['refl','vel', 'vel_texture_smoothed']
#  r_mom = ['refl','vel', 'vort_smooth', 'vel_texture']
#  r_mom = ['refl','vel', 'vort_smooth']
#  r_mom = ['refl','vel', 'vort', 'vort_smooth']
#  r_mom = ['refl','vel_smooth', 'zoom']
r_mom = ['refl','vel', 'Meso']
#  r_mom = ['refl','vel_smooth']
#  r_mom = ['refl','vel', 'vel_savgol_axis0']
#  r_mom = ['refl','vel','vel_savgol', 'vel_savgol_axis1']

#list; the timeseries to plot (can include multiple) (leave list empty to plot only radar) 
    #  histogram, wind, Thetae, Thetav, tfast
#  Line_Plot= []
#  Line_Plot= ['clicker']
#  Line_Plot= ['histogram']
#  Line_Plot= ['wind','Thetav']
Line_Plot= ['Thetav']

# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## What type of radar do you want to plot 
#  Radar_Plot_Type = 'NOXP_Plotting'
Radar_Plot_Type = 'KA_Plotting'
#  Radar_Plot_Type = 'WSR_Plotting'


# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## Timeframe
# ***********
## List of days to plot
#  day_list= ['20190608']
#  day_list= ['20190523']
day_list= ['20190524']
#  day_list= ['20190517','20190518','20190520']# '20190523', '20190524','20190525','20190526',
#  day_list= ['20190520', '20190523', '20190524','20190525','20190526',
#  day_list= ['20190517','20190518','20190520', '20190523', '20190524','20190525','20190526',
             #  '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190518','20190520', '20190523', '20190524','20190525','20190526',
             #  '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190525','20190526','20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190523', '20190524','20190525','20190526','20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190608', '20190611', '20190613', '20190615']

## Crop the start or end time to a time you specify (comment out to use datafrom the whole day)
#  tstart = dt.datetime(int(day_list[0][0:4]), int(day_list[0][4:6]), int(day_list[0][6:8])+1, 1, 32, 0)
#  tend = dt.datetime(int(day_list[0][0:4]), int(day_list[0][4:6]), int(day_list[0][6:8]), 23, 40, 0)


# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## Data controls  
# *******************************************
#whether to bother trying to read in this datatype at all (timesaver) 
Read_Data = {'Meso': True, 'Radar':True}


# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## Plot visual controls (aka what is plotted)
# *******************************************
#Time series controls 
    # X: restrict the time series to display only a timespan +/- 2x time from the scantime 
    # Y: hardcode a y axis extent for the timeseries
    # Wind_Pform: if you are plotting a wind time series what platform would you like it to be 
    # Mask: include in the list the type names for the pforms you would like masked 
lineplt_control = { 'Axis_Extents':{'X':30, 'Y':None}, #either x or y can be set to None for full range
                    'Wind_Pform':'Prb2',
                    #  'Wind_Pform':'CoMeT3',
                    'Mask': ['NSSL'],
                    'Deriv': True,
                    #  'H':{'var':'vort', 'xmin':-1, 'xmax':1}
                    'H':{'var':'vel_smooth', 'xmin':False, 'xmax':False}
                    }

Surge_controls= { 'ray_selection': 10, #degrees off of 90
                  'offset_dist': [500, 1000, 2000], #dist in m
                  'Make_new_pnts': False,
                  'Read_surge_stuff':True 
                 }
#Radar plot controls
    # p_tilt: what radar elevation tilt (deg) do you want to plot
    # Centered_Pform: what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    # offsetkm: size of the plot domain (times by 2): 21 is the best for KA
radar_controls= {'layout':{"Rows":1, "Cols":None},
                 #  'p_tilt': [0.5],
                 #  'p_tilt': [0.5, 1.0,1.5],
                 'p_tilt': [1.0],#, 1.5],
                 #  'Centered_Pform': 'WinS',
                 #  'Centered_Pform': 'CoMeT2',
                 #  'Centered_Pform': 'Prb1',
                 #  'Centered_Pform': 'Ka2',
                 'Centered_Pform': 'P_Radar',
                 'offsetkm':{'KA_Plotting':21, 'NOXP_Plotting':20, 'WSR_Plotting':40, 'Zoom':'Fitted'},
                 #  'offsetkm':{'KA_Plotting':15, 'NOXP_Plotting':10, 'WSR_Plotting':40},
                 'distance_indicator':{'square_grid':True, 'range_rings':True, 'ring_int':1000}, 
                 'colorbars':{'thermo':True, 'mom':True},
                 'smooth':13
                 }

#what to plot on top of the radar image
    # Marker: the marker on the radar plot 
    # Label: a text label on the radar plot 
    # Colorline: the cline that extends +- x min from a platfrom maker (also controls width of greybox on tseries)
    # Colorline: var: which var to plot on the colorline overlay (current options; Thetae, Thetav)
overlays = {'KA': {'Marker': True, 'Label': False, 'RHI_ring': True, 'RHI_Lab': False},
            'WSR': {'Marker': True, 'Label': False},
            'NOXP': {'Marker': True, 'Label': False},
            'MAINR':{'Marker':False, 'Label': False}, # always leave this one False
            'NSSL': {'Marker': True, 'Label': False, 'qcflags':['qc1', 'qc2', 'qc3', 'qc4']},
            'UNL': {'Marker': True, 'Label': False},
            'STN_I': {'Marker': True, 'Label': False},
            'UAS': {'Marker': False, 'Label': False},
            'Colorline': {'NSSL':True, 'UNL':True, 'UAS':False, 'cline_extent':5, 'Var':"Thetav", 'Pert':True},
            'LSR': {'Marker': False, 'Label': False},
            'Tracks': {'Marker': False, 'Label': False},
            'Contour':{'Lines': False, 'Interval':5, 'Label':True, 'Var':'vel_fix'},
            'surge_lines':{'Marker': True, 'Lines':True, 'zoom_boxes':True, 'zoom_only':True},
            'GEO': { 'OX': False, 'country_roads': False, 'hwys': False, 'county_lines': False, 'state_lines': False}
           }


# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## Troubleshooting (True/False variable)
# *****************
#  Set to True to print statements when entering/leaving definitions
    ## (helpful to understanding whats going on in the code but prints alot of lines to your terminal)
print_long = False 
#  Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting
    ## there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = False
#  Print radar metadata
print_radar_info=False 


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
