#import needed modules
######################
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patheffects as PathEffects
from matplotlib.ticker import (LinearLocator, FixedLocator, MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as dt
from datetime import datetime, date, timedelta
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import ShapelyFeature,NaturalEarthFeature
from cartopy.io.shapereader import Reader
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
import metpy
from metpy.plots import ctables
from metpy.units import units
import metpy.calc as mpcalc
import pandas as pd
import numpy as np
import numpy.ma as ma
import xarray as xr
import osmnx as ox
import os, os.path
from scipy import ndimage, interpolate
from operator import attrgetter
from collections import namedtuple
import pyart, sys, traceback, shutil, glob, gc
import time
import nexradaws
from pathlib import Path
from os.path import expanduser
from joblib import Memory


## Imports form other files
############################
import config #this is the file with the plotting controls to access any of the vars in that file use config.var

#rename a few commonly used vars so that the config.var does not have to be used repeatedly
print_long, e_test = config.print_long, config.e_test

## Read in defns that I have stored in another file (for ease of use/consistancy accross multiple scripts)
from shared_defns import Add_to_DATA, pform_names, error_printing, Platform, Radar, R_Plt, Torus_Insitu, Stationary_Insitu

##########################################################################
###########################
## PLOTTING DEFINTIONS  ###
###########################

def ppiplot(Data, print_long, e_test):
    ''' Initial plotting defenition: sets up fig size, layout, font size etc and will call timeseries and radar subplots
    ----------
    INPUTS
    ka_info : dictionary
        contains info about the ka radars (ie time, locations, sweep, azimuth, and RHI angles)
    var: dictionary
        contains info about the variable being overlayed (ie Thetae etc). Contains name and the max and min value
    pform: dictionary
        contains the pandas dataframes of the initu platforms
    #t_plotting, r_plotting : string (description not right)
    #    true/false variable that indicates whether only radar should be plotted (aka no Timeseries if True)
    #    (I have not actually tried this yet .... but in theroy this should work (you would need to play with formating))
    filesys & day: string
        base path to data and strorage, day of deployment YYYYMMDD
    print_long & e_test: strings
        True/False vars that control how much information is printed out to the terminal window. Helpful for debugging but can be overkill
    '''
    if print_long == True: print('~~~~~~~~~~~Made it into ppiplot~~~~~~~~~~~~~~~~~~~~~')
    #  print('******')
    #  print(matplotlib.rcParams)
    #  print('******')
    SMALL_SIZE, MS_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 23, 28, 33, 50
    plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=MS_SIZE)       # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    ## Set up plot size
    if config.t_plotting == False:
        fig = plt.figure(figsize=(20, 10), facecolor='white')
        gs = GridSpec(nrows=2, ncols=1)
    else:
        ## Set up figure
        fig = plt.figure(figsize=(32, 20), facecolor='white')
        ## Establish the gridspec layout
        gs = GridSpec(nrows=2, ncols=5, width_ratios=[.5, 8, 1, 8, .25], height_ratios=[3, 2], wspace=.25, hspace=.1)
        ## There are extra columns which are spacers to allow the formating to work
        # Uncomment to visualize the spacers
        #  ax_c = fig.add_subplot(gs[0, 2])
        #  ax_l = fig.add_subplot(gs[0, 0])
        #  ax_s = fig.add_subplot(gs[0, 4])

    ## Plot each radar subplots #only set up for two radar subplots currently
    t_R = R_Plt(Data)
    col = -1
    for mom in config.r_mom[:]:
        # row and col pos of the subplot in the gridspec
          # row should always be topmost row (row= 0)
          # first (left) subplot should be col=1 & the second (right) subplot should be col=3 due to the extra spacing "subplots" in the gridpec
        row, col = 0, col+2
        print('Radar moment: '+mom+', Gridspec row: '+str(row)+', GridSpec col: '+str(col))
        #if you are plotting the left subplot include the legend
        if col == 1: leg= True
        else: leg= False

        ## Make the subplot
        #### * * * * * * * *
        post = radar_subplots(mom, Data, t_R, fig, gs[row, col], leg, print_long, e_test)

    ## Make the Times Series Subplots (if applicable)
    #### * * * * * * * * * * * * * * * * * ** * * * *
    if config.t_plotting == True: time_series(Data, fig, gs, print_long, e_test)

    ## Plot platform colorbar
    cbar_ax = plt.axes([.514, post.y0, .014, post.y1-post.y0])#left, bottom, width, height
    cbar = plt.colorbar(Data['Var'].CS3, cax=cbar_ax, orientation='vertical', label=Data['Var'].v_lab, ticks=MaxNLocator(integer=True))#,ticks=np.arange(Data['p_var'].global_min, Data['p_var'].global_max+1,2))

    ## Plot title
    plt.suptitle(Data['P_Radar'].name+' '+str(config.p_tilt)+r'$^{\circ}$ PPI '+Data['P_Radar'].fancy_date_str, y=.92)

    ## Finish plot
    if config.wind == True: 
        output_name = config.filesys+'TORUS_Data/'+config.day+'/mesonets/plots/'+Data['P_Radar'].name+'_wind_'+config.p_var+'_'+Platform.Scan_time.strftime('%m%d_%H%M')+'.png'
    else: 
        output_name = config.filesys+'TORUS_Data/'+config.day+'/mesonets/plots/'+Data['P_Radar'].name+'_'+config.p_var+'_'+Platform.Scan_time.strftime('%m%d_%H%M')+'.png'
    print(output_name)
    start_time = time.time()
    plt.savefig(output_name, bbox_inches='tight', pad_inches=.3)
    print("My line took", time.time() - start_time, "to run")
    plt.close()

    ## Makes a ding noise
    print('\a')
    if print_long == True: print('~~~~~~~~~~~made it through ppiplot~~~~~~~~~~~~~~~~~~')
    print('Done Plotting \n ***************************************************************************************************')

# * * * * * *  *
def radar_subplots(mom, Data, t_R, fig, sub_pos, leg, print_long, e_test):
    ''' Plots each of the radar subplots and calls for the markers to be plotted for all the additional platforms
    ----
    INPUTS: mom: str, indicates which radar moment you would like to plot (ie Reflectivity, Velocity, Spectrum Width etc)
            fig, gs, & display: information about the plot
            ka_info, var, & pform: dictionarys, as described in the ppiplot comment
    '''
    if print_long == True: print('~~~~~~~~~~~made it into radar_subplots~~~~~~~~~~~~~~')

    ## SET UP VARS FOR EACH RADAR MOMENTS
    if mom == 'refl':
        p_title, field, c_label, c_scale, vminb, vmaxb = 'Reflectivity', 'refl_fix', 'Radar Reflectivity [dbz]', 'pyart_HomeyerRainbow', -30., 30.
    elif mom == 'vel':
        p_title, field, c_label, c_scale, vminb, vmaxb = 'Radial Velocity', 'vel_fix', 'Velocity [m/s]', 'pyart_balance', -40., 40.
    else: print("Hey what just happened!\n Check the Radar moments for spelling")


    ## SET UP SUBPLOTS
    ax_n = fig.add_subplot(sub_pos, projection= t_R.R_Proj)
    ax_n.text(.5, -.065, p_title, transform= ax_n.transAxes, horizontalalignment='center', fontsize=40) #the radar subplot titles
    t_R.display.plot_ppi_map(field, Data['P_Radar'].swp, title_flag=False, colorbar_flag=False, cmap=c_scale, ax=ax_n, vmin=vminb, vmax=vmaxb, min_lat=t_R.Domain.ymin, max_lat=t_R.Domain.ymax, min_lon=t_R.Domain.xmin, max_lon=t_R.Domain.xmax, embelish=False)

    ## PLOT PLATFORMS AS OVERLAYS(ie marker,colorline etc) ON RADAR
    #  iterate over each object contained in dict Data (returns the actual objects not their keys)
    legend_elements = [] #empty list to append the legend entries to for each platfrom that is actually plotted
    for p in Data.values():
        #  print(vars(p))
        ######
        #Plot Inistu Torus Platforms (if desired and available)
        if isinstance(p, Torus_Insitu):
            # * * *
            if config.NSSLm == True:
                if p.type == 'NSSL':
                    if print_long == True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    p.plot_Tpform(Data, ax_n, print_long, e_test)
            # * * *
            if config.NEBm == True:
                if p.type == 'UNL':
                    if print_long == True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    p.plot_Tpform(Data, ax_n, print_long, e_test)
            # * * *
            if config.UASm == True:
                if p.type == 'UAS':
                    if print_long == True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    p.plot_Tpform(Data, ax_n, print_long, e_test, border_c='xkcd:very pale green', labelbias=(0,0.01))

        ######
        #Plot Stationary Inistu Platforms (if desired and available)
        if isinstance(p, Stationary_Insitu):
            # * * *
            if config.MESONETSm == True:
                if p.type == 'WTM':
                    if print_long == True: print(p.name)
                    # determine if any of the sites fall within the plotting domain
                    sites_subdf, valid_sites = p.grab_pform_subset(print_long, e_test, Data, bounding= t_R.Domain)
                    # if there are sites within the domain plot the markers and include in the legend
                    if valid_sites == True:
                        legend_elements.append(p.leg_entry)
                        ax_n.plot(sites_subdf.lon, sites_subdf.lat, transform=ccrs.PlateCarree(), marker=p.m_style, linestyle='None', markersize=23, color=p.m_color)
                        # include labels for the sites on the plots
                        if config.WTxM_lab == True:
                            for x, y, lab in zip(sites_subdf['lon'], sites_subdf['lat'], sites_subdf['Stn_ID']):
                                trans = t_R.R_Proj.transform_point(x+.009, y-.002, ccrs.Geodetic())
                                ax_n.text(trans[0], trans[1], lab, fontsize= 20)
            # * * *
            if config.METARm == True:
                if p.type == 'METAR':
                    if print_long == True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    print("To Be filled in ")
            # * * *
            if config.ASOS_AWOSm == True:
                if p.type == 'ASOS/AWOS':
                    if print_long == True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    print("To Be filled in ")

        #####
        #Plot Radar Platforms (if desired and available)
        if isinstance(p, Radar):
            # * * *
            if config.KAm == True:
                if p.type == 'KA':
                    if print_long == True: print(p.name)
                    #det if the radar is located within the area included in the plot at the time of the plot
                    p_deploy = p.grab_pform_subset(print_long, e_test, Data, bounding= t_R.Domain, Single_Point= True)

                    #if these conditions are met then plot the radar marker
                    if p_deploy == True:
                        legend_elements.append(p.leg_entry)
                        ## Plot the marker
                        ax_n.plot(p.lon, p.lat, marker=p.m_style, transform=ccrs.PlateCarree(), color=p.m_color, markersize=18,
                                  markeredgewidth=5, path_effects=[PathEffects.withStroke(linewidth=15, foreground='k')], zorder=10)
                        ## Plot RHI spokes
                        if config.rhi_ring == True: p.rhi_spokes_rings()
                        ## Optional textlabel on plot
                        if config.KA_lab == True:
                            ax_n.text(p.lon+.009, p.lat-.002, p.name, transform=ccrs.PlateCarree(), path_effects=[PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')])
            # * * *
            if config.WSRm == True:
                if p.type == 'WSR':
                    if print_long == True: print(p.name)
                    # determine if any of the sites fall within the plotting domain
                    sites_subdf, valid_sites = p.grab_pform_subset(print_long, e_test, Data, bounding= t_R.Domain)
                    # if there are sites within the domain plot the markers and include in the legend
                    if valid_sites == True:
                        legend_elements.append(p.leg_entry)
                        ax_n.plot(sites_subdf.lon, sites_subdf.lat, transform=ccrs.PlateCarree(), marker=p.m_style, linestyle='None', markersize=20, color=p.m_color)
                        # include labels for the sites on the plots
                        if config.WSR88D_lab == True:
                            for x, y, lab in zip(sites_subdf['lon'], sites_subdf['lat'], sites_subdf['Stn_ID']):
                                trans = t_R.R_Proj.transform_point(x, y, ccrs.Geodetic())
                                ax_n.text(trans[0]+.02, trans[1]-.01, lab, fontsize= 27)
            # * * *
            if config.NOXPm == True:
                if p.type == 'NOXP':
                    if print_long == True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    print('Code not written yet')

    ## DEAL WITH COLORBARS
    # Attach colorbar to each subplot
    divider = make_axes_locatable(plt.gca())
    c_ax = divider.append_axes("right", "5%", pad="2%", axes_class=plt.Axes)
    sm = plt.cm.ScalarMappable(cmap=c_scale, norm=matplotlib.colors.Normalize(vmin=vminb, vmax=vmaxb))
    sm._A = []
    cb = plt.colorbar(sm, cax=c_ax, label=c_label)
    if leg == False:  #this means you are currently making the right subplot
        post = ax_n.get_position() ## Get position information to pass along to the remaining colorbar

    ## SET UP LEGENDS
    if leg == True: #this means you are currently making the left subplot
        #add legend for platform markers
        l = ax_n.legend(handles=legend_elements, loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(0,.5), handlelength=.1, title="Platforms", shadow=True, fancybox=True, ncol=1, edgecolor='black')
        l.get_title().set_fontweight('bold')
        post = np.nan

    ## PLOT BACKGROUND FEATURES
    if config.country_roads == True:
        ox.config(log_file=True, log_console=True, use_cache=True) #the config in this line has nothing to do with config.py
        G = ox.graph_from_bbox(t_R.Domain.ymax, t_R.Domain.ymin, t_R.Domain.xmax, t_R.Domain.xmin)
        ox.save_load.save_graph_shapefile(G, filename='tmp'+str(0), folder=config.filesys+'radarproc/roads/', encoding='utf-8')
        fname = config.filesys+'radarproc/roads/tmp'+str(0)+'/edges/edges.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(), ccrs.PlateCarree(), edgecolor='gray', linewidth=0.5)
        ax_n.add_feature(shape_feature, facecolor='none')
        shutil.rmtree(config.filesys+'radarproc/roads/tmp'+str(0)+'/')
    if config.hwys == True:
        fname = config.filesys+'radarproc/roads/GPhighways.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(), ccrs.PlateCarree(), edgecolor='grey')#edgecolor='black')
        ax_n.add_feature(shape_feature, facecolor='none')
    if config.county_lines == True:
        fname = config.filesys+'radarproc/roads/cb_2017_us_county_5m.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(), ccrs.PlateCarree(), edgecolor='gray')
        ax_n.add_feature(shape_feature, facecolor='none', linewidth=1.5, linestyle="--")
    if config.state_lines == True:
        states_provinces = cartopy.feature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='10m', facecolor='none')
        ax_n.add_feature(states_provinces, edgecolor='black', linewidth=2)
    if print_long == True: print('~~~~~~~~~~~Made it through radar_subplots~~~~~~~~~~~')

    return post

# * * * * * * *
def time_series(Data, fig, gs, print_long, e_test):
    ''' Plot the time series of p_var from the various instrument platforms
    ----
    INPUTS; fig & gs: relating to the plot layout
            ka_info, var, & pform: dictionarys, as described in the ppiplot comment
            r_plotting: True/False, If False will only plot timeseries (no radar).. In theroy have not actually done this yet
    '''
    if print_long == True: print('~~~~~~~~~~~Made it into time_series~~~~~~~~~~~~~~~~~')

    if config.r_plotting == True: ax_n= fig.add_subplot(gs[1, :])
    else: fig, ax_n = plt.subplots() #will only plot the time series (not the radar subplots)

    ## MAKE THE TIMESERIES
    #### * * * * * * * * *
    if config.wind == False:
        for p in Data.values():
            if isinstance(p, Torus_Insitu):
                if print_long == True: print('Plotting '+str(p.name)+' on time series')
                ## If the platform matches a type listed in TS_masking use masked data for the time series; else use the unmasked data
                if p.type in config.TS_masking[:]: plotting_data= p.ts_mask_df
                else: plotting_data= p.df[config.p_var].values

                ## Plot
                ax_n.plot(p.df['datetime'], plotting_data, linewidth=3, color=p.l_color, label=p.leg_str) #assigning label= is what allows the legend to work

        ## Set up XY axes tick locations
        ax_n.yaxis.set_major_locator(MultipleLocator(5)) # set up major tick marks (this is set up to go by 5's will want to change for diff vars)
        ax_n.yaxis.set_minor_locator(AutoMinorLocator(5)) # set up minor ticks (this have it increment by 1's will want to change for diff vars)

        ## Set up axes formats (lables/ ticks etc)
        ax_n.set_ylabel(Data['Var'].v_lab)

        ## Set up grid for plot
        ax_n.grid(which='major', axis='y', linestyle='--', linewidth=2)
        ax_n.grid(which='minor', axis='y', linestyle=':')

        ## Include the legend
        leg = ax_n.legend()
        for line in leg.get_lines(): line.set_linewidth(12)

    if config.wind == True:
        p = Data['Prb1']
        if print_long == True: print('Plotting '+str(p.name)+' on time series')
        l1 = ax_n.plot(p.df['datetime'], p.df['spd'], label='Wind Speed')
        ax_n.fill_between(p.df['datetime'], p.df['spd'], 0)

        ax2 = ax_n.twinx()
        l3 = ax2.plot(p.df['datetime'], p.df['dir'], '.k', linewidth=.5, label='Wind Direction')
        

        ax_n.set_ylabel('Wind Speed', multialignment='center')
        ax2.set_ylabel('Wind Direction (degrees)', multialignment='center')
        ax2.set_ylim(0, 360)
        #  ax2.set_yticks(np.arange(45, 405, 90), ['NE', 'SE', 'SW', 'NW'])
        ax_n.yaxis.set_major_locator(LinearLocator(numticks=9))
        # set up minor ticks 
        ax2.yaxis.set_minor_locator(FixedLocator(np.arange(45, 405, 90))) 
        # set up major tick marks 
        ax2.yaxis.set_major_locator(FixedLocator(np.arange(90, 360, 90)))
        ax2.grid(b=True, which='minor', axis='y', color='k', linestyle='--', linewidth=0.5)
        ax2.grid(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=1)
        ax2.tick_params(which='major', width=2, length=20, color='black')
        ax2.tick_params(which='minor', width=2, length=10, color='grey')

        lines = l1 + l3  
        labs = [line.get_label() for line in lines]
        #  ax2.legend(lines, labs, prop={'size': 12})
        ax2.legend(lines, labs)

    # if desired this subsets the timerange that is displayed in the timeseries 
    if config.ts_extent != None:
        ax_n.set_xlim(Platform.Scan_time - timedelta(minutes=config.ts_extent), Platform.Scan_time + timedelta(minutes=config.ts_extent))

    ax_n.margins(x=0, y=0)
    ax_n.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax_n.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M')) #strip the date from the datetime object
    ax_n.xaxis.set_minor_locator(AutoMinorLocator(6)) # set up minor ticks (should be a multiple of ten intervals (ie 10,20,30... min spans)
    ax_n.tick_params(which='major', width=2, length=20, color='black')
    ax_n.tick_params(which='minor', width=2, length=10, color='grey')
    ax_n.grid(which='major', axis='x', color='grey', linewidth=2.5)
    ax_n.grid(which='minor', axis='x')
    ax_n.set_xlabel('Time (UTC)')
    
    ## If makeing the timeseries in conjunction with radar subplots set up vertical lines that indicate the time of
    #  radar scan and the timerange ploted (via filled colorline) on the radar plots
    if config.r_plotting == True:
        ax_n.axvline(Platform.Scan_time, color='r', linewidth=4, alpha=.5, zorder=10)
        ax_n.axvspan(Platform.Scan_time - timedelta(minutes=config.cline_extent), Platform.Scan_time + timedelta(minutes=config.cline_extent), facecolor='0.5', alpha=0.4, zorder=10)

    if print_long == True: print('~~~~~~~~~~~Made it through time_series~~~~~~~~~~~~~~')

##############################################################################################
#######################
## RADAR DEFINITIONS ##
#######################


def get_WSR_from_AWS(start, end, radar_id, download_directory):
    '''
    Retrieve the NEXRAD files that fall within a timerange for a specified radar site from the AWS server
    ----------
    INPUTS radar_id : string;  four letter radar designation
           start: datetime;  start of the desired timerange
           end: datetime; end of the desired timerange
           temploc: string; location for the downloaded radarfiles
    -------
    RETURN radar_list : Py-ART Radar Objects
    '''

    #Create this at the point of use otherwise it saves everything and eventually crashes 
    conn = nexradaws.NexradAwsInterface()

    #Determine the radar scans that fall within the time randge for a given radar site
    scans = conn.get_avail_scans_in_range(start, end, radar_id)
    print("There are {} scans available between {} and {}\n".format(len(scans), start, end))
  
    #Dont download files that you alrady have ....
    path = config.filesys+ 'TORUS_Data/'+ config.day + '/radar/Nexrad/Nexrad_files'

    if not os.path.exists(path): 
        #  Path(path).mkdir(parents=True, exist_ok=True)
        Path(path).mkdir(parents=True)

    # missing_scans is a list of scans we don't have and need to download create_filepath returns tuple of 
    # (directory, directory+filename) [-1] returns the directory+filename
    missing_scans = list(filter(lambda x: not Path(x.create_filepath(path,False)[-1]).exists(), scans))

    # missing files is the list of filenames of files we need to download create_filepath returns tuple of 
    # (directory, directory+filename) [-1] returns the directory+filename
    missing_files_after = list(filter(lambda x: not Path(x.create_filepath(path,False)[-1]).exists(), scans))

    if len(missing_files_after) >0:
        print("ERROR: Some Radar Scans are Missing \n", missing_files_after)
        exit()

    radar_files = list(map(lambda x: x.create_filepath(path,False)[-1], scans))
    # Return list of files 
    return radar_files 


##############################################################################################
#####################################################
# Read in Data that does not change for each image ##
#####################################################

print('\n Read in Torus Platforms \n')
#print(pform_names('ALL')) #list of ALL possible platforms
subset_pnames = [] #array of platforms that actually have data for the time range we are plotting
Data = {} #dictionary in which all the class objects will be stored (will contain platform data, locations, plotting variables etc)

## Read in the data for the TORUS Insitu platforms (if available)
Data, subset_pnames = Add_to_DATA('TInsitu', Data, subset_pnames, print_long)
## Establish info for the plotting variable (aka max min etc)
Data, subset_pnames = Add_to_DATA('PVar', Data, subset_pnames, print_long)

print('\n Read in Stationary Platforms Arrays\n')
## Read in the data for the Stationary Array platforms (if available)
Data, subset_pnames = Add_to_DATA('STN_I', Data, subset_pnames, print_long)


###################################
# Create Plot for each radarfile ##
###################################
## If Radar will be plotted
if config.r_plotting == True:
    print('\n Yes Plot Radar \n')
    if config.KA_Plot == True:
        ## Get radar files
        fil = sorted(glob.glob(config.filesys+'TORUS_Data/'+config.day+'/radar/TTUKa/netcdf/*/dealiased_*'))

        ## Read in and check each radarfile in fil
        #### + + + + + + + + + + + + + + + + + + +
        for thefile in fil[:]:
            radar = pyart.io.read(thefile) ## Read the radar file

            ## We are currently not interested in plotting the RHI scans
            if radar.scan_type == 'rhi': print(str(thefile) + "\n This scan is an RHI \n **********************************")

            ## If file contains a ppi scan proceed
            elif radar.scan_type == 'ppi':
                n_swps = radar.nsweeps
                for swp_id in range(n_swps):
                    ## Det the actual tilt angle of a given sweep (returns an array)
                    tilt_ang = radar.get_elevation(swp_id)

                    ## Check to see if the radarfile matches the elevation tilt we are interested in
                    if np.around(tilt_ang[0], decimals=1) == config.p_tilt:
                        print(str(thefile) + "\n Producing Radar Plot:")

                        ## Read in radar data and add to Data dict
                        ##### + + + + + + + + + + + + + + + + + + +
                        #  Establish info for the main plotting radar (aka scantime etc) & and locations for other radars (if deployed)
                        Data, subset_pnames = Add_to_DATA('RADAR', Data, subset_pnames, print_long, MR_file=radar, swp=swp_id)
                        if print_long == True: print(Data, '\n')

                        ## Proceed to plot the radar
                        ##### + + + + + + + + + + + +
                        ppiplot(Data, print_long, e_test)

                    ## If the scans elevation tilt does not match the one we are plotting
                    else: print(str(thefile) + "\n Scan does not match the tilt currently being plotted \n Currently plotting: Elevation Tilt= "
                                 + str(config.p_tilt) +"\n This scan: Elevation tilt= "+ str(tilt_ang[0])+'\n**********************************')

    if config.WSR_Plot == True: 
        '''
        locate the nearest WSR88D site to the specified insitu instruments
        '''
        #find the locations of all WSR88D sites (outputs dict in format {Site_ID:{lat:...,lon:...,elav:...], ...})
        all_WSR = pyart.io.nexrad_common.NEXRAD_LOCATIONS
        #  print(json.dumps(all_WSR_locs, sort_keys=True, indent=4))

        #set up empty dataframe with site Ids as column names
        d_from_all_r = pd.DataFrame(columns = all_WSR.keys())

        #fill in dataframe with the distance from all 88D sites from each probe measurement
        for key in all_WSR:
            d_from_r = np.square(Data[config.Centered_Pform].df['lat']-all_WSR[key]['lat']) +np.square(Data[config.Centered_Pform].df['lon']-all_WSR[key]['lon'])
            d_from_all_r[key] = d_from_r

        #Determine which WS88D site is closest to the probe and add to the original probe dataframe
        Data[config.Centered_Pform].df['Radar_ID'] = d_from_all_r.idxmin(axis = 1)
        
        #Determine the unique radar sites to be plotted 
        r_ofintrest = Data[config.Centered_Pform].df.Radar_ID.unique()
        print(r_ofintrest)

        timeranges_each_r = pd.DataFrame()
        for r in r_ofintrest:
            print(r)
            t = Data[config.Centered_Pform].df.loc[Data[config.Centered_Pform].df.Radar_ID == r, ['datetime']].rename(columns={'datetime':r})
            ts, te = t.min(), t.max()
            timeranges_each_r = pd.concat([timeranges_each_r, t], axis=1)
            print("start ", ts[r])
            print("end ", te[r])
            print("***")

            radar_files = get_WSR_from_AWS(ts[r], te[r], r, config.filesys)

            print('********')
            print("Radar files to process:")
            print(radar_files)

            #open the downloaded files as pyart objects
            radar_list=[]
            i=0
            for radar_file in radar_files:
                i = i+1
                print("[{}] open_pyart, scan file_name = {}\n".format(i, radar_file))
                
                try: radar = cached_radar_from_nextrad_file (radar_file)
                except: print("Failed to convert file: ", radar_file)

                print(i,radar.info())

                #why am i calling radar and radarlist to this function when I am in a loop that is interating through these to values?
                ppiplot(r_only, radar, config.filesys, day, globalamin,globalamax, p_var, p_of_int,CS3, e_test)
        print(timeranges_each_r)


################################
# Create Timeseries only plot ##
################################
#Only plot timeseries (this code isn't fully fleshed out but in theroy this code is built in such a way to allow for this)
if config.r_plotting == False and config.t_plotting == True:
    print("Plot Timeseries only \n"+ str(config.filesys+'TORUS_Data/'+config.day+'/mesonets/NSSL/*.nc'))
    time_series(Data)
    fig.savefig('test2.png')
    plt.close()

###########################
print("ALL FINISHED")
