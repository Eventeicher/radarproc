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
p_azimuth= 1.0 #what radar azimuth scan do you want to plot
pri_tseries='p_var' #what would you like to plot in the first time series (if num_of_tseries=1 this is the only tseries to be plotted)
sec_tseries='wind' #what would you like in the second time series (only used if num_of_tseries =2)
# would you like to restrict the time series to display only a timespan +/- x time from the scantime (functionalitly not coded yet)
tseries_bound= False # to display full timespan set this var to False otherwise set this var to x (in min)
cline_extent= "placeholder" #how long would you like the colorlines for each platforms to cover (also control width of greybox on tseries) 

## Timeframe
# ***********
day = '20190524' #'YYYYMMDD'
#crop the start or end time to a time you specify (set to None to do the whole day)
tstart = dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),15,0,0)
tend = None

## Troubleshooting 
# *****************
#  Set to True to print statements when entering/leaving definitions
#  (helpful to understanding whats going on in the code but prints alot of lines to your terminal)
print_long = False #True/False variable
#  Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting
    ## there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = False #True/False variable

