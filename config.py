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
nCPU = -3
#  nCPU = 1
actually_plt_radar=True

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
#  r_mom = ['refl','vel', 'vel_texture', 'vel_test']
#  r_mom = ['refl','vel', 'vel_test']
#  r_mom = ['refl','vel', 'vel_texture','vel_texture_smoothed']
#  r_mom = ['refl','vel', 'vort_smooth', 'vel_texture']
#  r_mom = ['refl','vel', 'vort_smooth', 'zoom']
#  r_mom = ['refl','vel', 'vort', 'vort_smooth']
#  r_mom = ['refl','vel_smooth', 'zoom', 'Meso']
#  r_mom = ['refl','vel', 'Meso_azi', 'zoom']
#  r_mom = ['refl','vel', 'Meso_azi']
#  r_mom = ['refl','vel', 'Meso']
#  r_mom = ['refl','vel']
#  r_mom = ['refl','vel', 'MRMS']
r_mom = ['refl','vel', 'zoom']
#  r_mom = ['refl','vel_smooth', 'zoom']
#  r_mom = ['refl','vel_smooth', 'vel_test']
#  r_mom = ['refl','vel_smooth', 'Meso', 'zoom']
#  r_mom = ['refl','vel_smooth', 'Meso_test', 'zoom']
#  r_mom = ['refl','vel', 'Meso_azi']#, 'zoom']
#  r_mom = ['refl','vel', 'vel_savgol_axis0']
#  r_mom = ['refl','vel','vel_savgol', 'vel_savgol_axis1']

#list; the timeseries to plot (can include multiple) (leave list empty to plot only radar) 
    #  histogram, wind, Thetae, Thetav, tfast
#  Line_Plot= []
Line_Plot= ['clicker']
#  Line_Plot= ['histogram']
#  Line_Plot= ['wind','Thetav']
#  Line_Plot= ['Thetae','Thetav']

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
#  day_list= ['20190518']
#  day_list= ['20190524','20190527']
#  day_list= ['20190518','20190520']# '20190523', '20190524','20190525','20190526',
#  day_list= [ '20190608']
#  day_list= [ '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= [ '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190611', '20190613', '20190615']
#  day_list= [ '20190611', '20190613', '20190615']
#  day_list= ['20190608', '20190611', '20190613', '20190615']
#  day_list= [ '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190527','20190524','20190608', '20190517','20190518','20190520', '20190523','20190525',
             #  '20190528', '20190611', '20190613', '20190615']
#  day_list= ['20190611', '20190613', '20190615']
#  day_list= ['20190611', '20190613', '20190615']
day_list= [ '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190517','20190518','20190520', '20190523', '20190524','20190525','20190526',
             #  '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']

## Crop the start or end time to a time you specify (comment out to use datafrom the whole day)
#  tstart = dt.datetime(int(day_list[0][0:4]), int(day_list[0][4:6]), int(day_list[0][6:8]), 23,11, 0)
#  tend = dt.datetime(int(day_list[0][0:4]), int(day_list[0][4:6]), int(day_list[0][6:8]), 1, 5, 0)


# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # #
## Data controls  
# *******************************************
#whether to bother trying to read in this datatype at all (timesaver) 
Read_Data = {'Meso': True, 'Radar':True}
Pert_times= { '20190517': {'Prb2':{'starth':21, 'startm':40, 'endh':22, 'endm':30}}, #maybe 35
              '20190518': {'All':{'starth':21, 'startm':30, 'endh':22, 'endm':20}},
              '20190520': None,
              '20190523': {},
              '20190524': {'All except p1':{'starth':18, 'startm':10, 'endh':18, 'endm':55}},
              '20190525': {'Prb1':{'starth':22, 'startm':35, 'endh':23, 'endm':10}},
              '20190526': None,
              '20190527': {'All':{'starth':20,'startm':0, 'endh':20, 'endm':15}},
              '20190528': None,
              '20190608': {'All':{'starth':22, 'startm':30, 'endh':23, 'endm':20}},
              '20190611': {'All': {'starth':23, 'startm':5, 'endh':0, 'endm':0 }},
              '20190613': {'All': {'starth':22, 'startm':45, 'endh':0, 'endm':30 }},
              '20190615': {'Prb1':{'starth':21, 'startm':35, 'endh':22, 'endm':10}}}

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
                    'Deriv': False, #thermo derivative
                    'Deriv_Rdata': True,
                    #  'H':{'var':'vort', 'xmin':-1, 'xmax':1}
                    'H':{'var':'vel_smooth', 'xmin':False, 'xmax':False},
                    'smooth_window':17, #17 None or odd number
                    'ray_style':'averagedrays' #allrays, pntrays, averagedrays ##to use while plotting surge data
                    }

Surge_controls= {'Surges': {'ray_selection_thresh': 20, #degrees off of 90
                            'ray_selection_method': 'local', # 'average' or 'local' referring to the type of slope intersection
                            'offset_dist': [2000] #dist in m
                            #  'offset_dist': [500, 1000, 2000] #dist in m
                            #  'offset_dist': [2000] #dist in m
                             },
                 'Mesos':  {'radius': 4023},
                 'Zoom': {'mom':'vel_smooth'},
                 'Feature_IDing': {'Make_new_pnts': {'activate_clicker':False, 'Type':'Meso'},#Surge'},   
                                   'Existing_pnts': {'Read_in_surges': True, 'Read_in_mesos':False,
                                                     'cross_rpforms':{'allow':False, 'within_time':5, 'tilt_allowance':.5}
                                                     }
                                    }
                 }
Conv_controls={'Valid_Rays': True, #should you subset to only the rays that fall in valid range when making the con plots 
              'add_to_con_csv':False,
              'make_scatterplts':False 
                }
#Radar plot controls
    # p_tilt: what radar elevation tilt (deg) do you want to plot
    # Centered_Pform: what should the radar subplots be centered on (for the plotting radar use 'P_Radar')
    # offsetkm: size of the plot domain (times by 2): 21 is the best for KA
radar_controls= {'layout':{"Rows":1, "Cols":None},
                 #  'p_tilt': [0.5, 1.0],
                 #  'p_tilt': [0.5, 1.0,1.5],
                 'p_tilt': [1.0],#, 1.5],
                 #  'Centered_Pform': 'WinS',
                 #  'Centered_Pform': 'CoMeT1',
                 #  'Centered_Pform': 'Ka2',
                 #  'Centered_Pform': 'Prb1',
                 'Centered_Pform': 'P_Radar',
                 'offsetkm':{'KA_Plotting':21, 'NOXP_Plotting':40, 'WSR_Plotting':30, 'Zoom':'Fitted'},
                 #  'offsetkm':{'KA_Plotting':15, 'NOXP_Plotting':10, 'WSR_Plotting':40},
                 'distance_indicator':{'square_grid':True, 'range_rings':True, 'ring_int':1000}, 
                 'colorbars':{'thermo':False, 'mom':True},
                 'smooth':13,
                 'Use_downloaded_files': True
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
            'Colorline': {'NSSL':False, 'UNL':False, 'UAS':False, 'cline_extent':5, 'Var':'Thetav', 'Pert':False},
            'LSR': {'Marker': False, 'Label': False},
            'Tracks': {'Marker': False, 'Marker_test':False, 'Label': False},
            'Contour':{'Lines': False, 'Interval':5, 'Label':True, 'Var':'vel_fix'},
            'surge_lines':{'Marker': True, 'Lines':True, 'zoom_boxes':True, 'zoom_only':True},
            'Mesos':{'Marker': True},
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
#g_root = '/Volumes/Samsung_T5/Research' # Ellie External SSD
g_root = '/media/chris/Samsung_T5/Research' # Ellie External SSD

g_TORUS_directory = g_root + '/TORUS_Data/'   # Our Data
g_csv_directory = g_root+ '/csv_files/'
g_cache_directory = g_root + '/cache/' # Cache
g_download_directory = g_cache_directory + 'downloads/' # Download Storage
g_roads_directory = g_root + '/roads/' # Roads data

# Setup Function Cache for speedup
if not os.path.exists(g_cache_directory): Path(g_cache_directory).mkdir(parents=True, exist_ok=True)
