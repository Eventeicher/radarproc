#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import needed modules
######################
import datetime as dt

# # # # # # # # # # # #  # # #  # # #  # # # # # #  # # # # # # # # # # # 
## VARIABLES ##
###############
## File paths
# ************
filesys='/Users/severe2/Research/'
temploc='/Volumes/Samsung_T5/Research/TORUS_Data/'
KA_Plot, WSR_Plot = False, True

## Timeframe
# ***********
day = '20190524' #'YYYYMMDD'
#crop the start or end time to a time you specify (set to None to do the whole day)
#  tstart = dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),20,15,0)
#  tend = dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),22,45,0)
#  ts_tend = dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),20,45,0)
tstart, tend= None, None

## Plot layout controls 
# ***********************
#  (some of these options will need to be expanded upon to work.... 
#  but hypothetically the code is set up to allow these options(with some modifications))
r_plotting= True #would you like to plot radar images? (set to False to only plot timeseries)
t_plotting= True #would you like to plot the timeseries? (set to False to only plot radar)
num_of_tseries= 1 #how many timeseries would you like included in layout (I envision a max of 2)

## Plot visual controls (aka what is plotted) 
# *******************************************
p_var = "Thetav" #which var to plot (current options; Thetae, Thetav)
offsetkm = 21
Centered_Pform = 'Prb1' #what should the radar subplots be centered on (for the plotting radar use 'P_Radar') 
p_tilt= 1.0 #what radar elevation tilt (deg) do you want to plot
rhi_ring= True #do you want the rhi spokes for the KA radars
r_mom =['refl','vel'] #list; the radar moments to plot on the left and right subplots respectively (current options are 'refl' or 'vel)

#which other platforms do you wish to plot as overlays to the radar image (aka as markers, colorlines etc) if data is available 
KAm, WSRm, NOXPm, NSSLm, NEBm, UASm, ASOS_AWOSm, MESONETSm, METARm,= True, True, False, True, True, False, False, True, False 
WTxM_lab, WSR88D_lab, KA_lab = True, True, False
cline_extent= 5 #how long would you like the colorlines for each platforms to cover +-x min (also control width of greybox on tseries) 
ts_extent= 30 #so actually 60 min 
country_roads, hwys, county_lines, state_lines = False, False, False, False #background plot features

#list; include in the list the type names for the pforms you would like masked (this only controls whats plotted on the Tseries) 
#  if none leave as blank list, will not work if a mask had not been previously applied to the dataset (which can be done in the classes)
TS_masking= ['NSSL']
wind = False
pri_tseries='p_var' #what would you like to plot in the first time series (if num_of_tseries=1 this is the only tseries to be plotted)
sec_tseries='wind' #what would you like in the second time series (only used if num_of_tseries =2)
# would you like to restrict the time series to display only a timespan +/- x time from the scantime (functionalitly not coded yet)
tseries_bound= False # to display full timespan set this var to False otherwise set this var to x (in min)



## Troubleshooting 
# *****************
#  Set to True to print statements when entering/leaving definitions
#  (helpful to understanding whats going on in the code but prints alot of lines to your terminal)
print_long = True #True/False variable
#  Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting
    ## there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = False #True/False variable

