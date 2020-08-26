#import needed modules
######################
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patheffects as PathEffects
from matplotlib.ticker import (LinearLocator, FixedLocator, MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib import ticker
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
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
from joblib import Memory, Parallel, delayed
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# To run with a) warnings and b) stack trace on abort
# python3 -Walways  -q -X faulthandler plot_nexrad_insitu.py

totalcompT_start = time.time()
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
def ppiplot(Data, print_long, e_test, start_comptime):
    ''' Initial plotting defenition: sets up fig size, layout, font size etc and will call timeseries and radar subplots
    ----------
    INPUTS: 
    Data: (dictionary) 
        contains objects corresponding to all available datasets (radar, torus insitue platforms etc) and  which particular variable 
        should be plotted on time series and colorlines 
    print_long & e_test: (bool strings) 
        True/False vars that control how much information is printed out to the terminal window. Helpful for debugging but can be overkill
    start_comptime: 
        Time at which you first begain plotting this particular image (will help to report out how long it took to create image)
    '''
    if print_long == True: print('~~~~~~~~~~~Made it into ppiplot~~~~~~~~~~~~~~~~~~~~~')
    SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 15, 23, 50
    plt.rc('font', size=SMALL_SIZE)         # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)       # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    ## Set up plot size
    if config.t_plotting == False:
        fig = plt.figure(figsize=(20, 10), facecolor='white')
        gs = GridSpec(nrows=2, ncols=1)
    else:
        ## Set up figure
        fig = plt.figure(figsize=(32, 20), facecolor='white')
        ## Establish the gridspec layout
        if len(config.Time_Series) == 1:
            gs = GridSpec(nrows=2, ncols=5, width_ratios=[.5, 8, 1, 8, .25], height_ratios=[3, 2], wspace=.25, hspace=.1)
            ts_gs = GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[1, :])
        if len(config.Time_Series) == 2:
            gs = GridSpec(nrows=2, ncols=4, width_ratios=[.15,8, 1, 8], height_ratios=[4,3], hspace=.1)
            ts_gs = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, :], height_ratios=[2, 3], hspace=0)
        ## There are extra columns which are spacers to allow the formating to work
        # Uncomment to visualize the spacers
        #  ax_c = fig.add_subplot(gs[0, 2])
        #only uncomment thes two if you are only doing 1 time series
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
        radar_subplots(mom, Data, t_R, fig, gs[row, col], leg, print_long, e_test)

    ## Make the Times Series Subplots (if applicable)
    #### * * * * * * * * * * * * * * * * * ** * * * *
    if config.t_plotting == True: 
        row = -1
        for ts in config.Time_Series[:]:
            row, col = row + 1, ':' 
            print('Time Series plot: '+ ts +', Gridspec row: '+str(row)+', GridSpec col: '+str(col))
            time_series(ts, Data, fig, ts_gs[row, :], print_long, e_test)
    #  ts_gs.align_ylabels()

    ## Plot title
    plt.suptitle(Data['P_Radar'].site_name+' '+str(config.p_tilt)+r'$^{\circ}$ PPI '+Data['P_Radar'].fancy_date_str, y=.92)
    file_string = '_'.join(config.Time_Series)

    ## Finish plot
    if Data['P_Radar'].name in pform_names('KA'):
        output_name = config.temploc+config.day+'/mesonets/plots/KA/'+Data['P_Radar'].site_name+'_'+Platform.Scan_time.strftime('%m%d_%H%M')+'_'+file_string+'.png'
    if Data['P_Radar'].name == 'WSR88D':
        output_name = config.temploc+config.day+'/mesonets/plots/WSR/'+Platform.Scan_time.strftime('%m%d_%H%M')+'_'+file_string+'_'+Data['P_Radar'].site_name+'.png'
    print(output_name)

    plt.savefig(output_name, bbox_inches='tight', pad_inches=.3)
    print("Plot took "+ str(time.time() - start_comptime)+ " to complete")
    plt.close()

    ## Makes a ding noise
    print('\a')
    if print_long == True: print('~~~~~~~~~~~made it through ppiplot~~~~~~~~~~~~~~~~~~')
    print('Done Plotting \n \n***************************************************************************************************')

# * * * * * *  *
def radar_subplots(mom, Data, t_R, fig, sub_pos, leg, print_long, e_test):
    ''' Plots each of the radar subplots including the marker overlays corresponding to the additional platforms
    ----
    INPUTS: 
    mom: str, indicates which radar moment you would like to plot (ie Reflectivity, Velocity, Spectrum Width etc)
    Data: dict as described in the ppiplot defn 
    t_R: .....fill in objected containing subplot info such as domain and radar.display
    fig, sub_pos: info about the positioning of the subplot within the overall layout of the figure (sub_pos is the gridspec pos)
    leg: bool str whether or not you want a legend associated with this particular subplot
    print_long & e_test: bool str as described in the ppi defn 
    '''
    if print_long == True: print('~~~~~~~~~~~made it into radar_subplots~~~~~~~~~~~~~~')

    ## SET UP VARS FOR EACH RADAR MOMENTS
    if mom == 'refl':
        p_title, c_label, c_scale = 'Reflectivity', 'Radar Reflectivity [dbz]', 'pyart_HomeyerRainbow'
        if Data['P_Radar'].name in pform_names('KA'):
            field, vminb, vmaxb, sweep = 'refl_fix', -30., 30., Data['P_Radar'].swp
        if Data['P_Radar'].name == 'WSR88D':
            field, vminb, vmaxb, sweep =  'reflectivity', -10., 75., Data['P_Radar'].swp[0]
    elif mom == 'vel':
        p_title, c_label, c_scale = 'Radial Velocity', 'Velocity [m/s]', 'pyart_balance'
        if Data['P_Radar'].name in pform_names('KA'):
            field, vminb, vmaxb, sweep = 'vel_fix', -40., 40., Data['P_Radar'].swp
        if Data['P_Radar'].name == 'WSR88D':
            field, vminb, vmaxb, sweep = 'velocity', -40., 40., Data['P_Radar'].swp[1]
    else: print("Hey what just happened!\n Check the Radar moments for spelling")


    ## SET UP SUBPLOTS
    ax_n = fig.add_subplot(sub_pos, projection= t_R.R_Proj)
    ax_n.text(.5, -.065, p_title, transform= ax_n.transAxes, horizontalalignment='center', fontsize=40) #the radar subplot titles
    try:
        t_R.display.plot_ppi_map(field, sweep, title_flag=False, colorbar_flag=False, cmap=c_scale, ax=ax_n, vmin=vminb, vmax=vmaxb, min_lat=t_R.Domain.ymin, max_lat=t_R.Domain.ymax, min_lon=t_R.Domain.xmin, max_lon=t_R.Domain.xmax, embelish=False)
    except: 
        print(Data['P_Radar'].rfile.info())
        print(Data['P_Radar'].Scan_time)
        #  exit()
    t_R.display.label_xaxis_x()

    ## PLOT PLATFORMS AS OVERLAYS(ie marker,colorline etc) ON RADAR
    #  iterate over each object contained in dict Data (returns the actual objects not their keys)
    legend_elements = [] #empty list to append the legend entries to for each platfrom that is actually plotted
    for p in Data.values():
        #  print(vars(p))
        
        ######
        #Plot Inistu Torus Platforms (if desired and available)
        if isinstance(p, Torus_Insitu):
            if print_long == True: print(p.name)
            legend_elements.append(p.leg_entry)
            p.plot_Tpform(Data, ax_n, print_long, e_test)
            #  p.plot_Tpform(Data, ax_n, print_long, e_test, border_c='xkcd:very pale green', labelbias=(0,0.01))

        ######
        #Plot Stationary Inistu Platforms ( aka mesonets and ASOS) if desired and available
        if isinstance(p, Stationary_Insitu):
            if print_long == True: print(p.name)
            # determine if any of the sites fall within the plotting domain
            sites_subdf, valid_sites = p.grab_pform_subset(print_long, e_test, Data, bounding= t_R.Domain)
            # if there are sites within the domain plot the markers and include in the legend
            if valid_sites == True:
                legend_elements.append(p.leg_entry)
                ax_n.plot(sites_subdf.lon, sites_subdf.lat, transform=ccrs.PlateCarree(), marker=p.m_style, linestyle='None', markersize= p.m_size, color=p.m_color)
                # include labels for the sites on the plots
                if p.marker_label == True:
                    for x, y, lab in zip(sites_subdf['lon'], sites_subdf['lat'], sites_subdf['Stn_ID']):
                        trans = t_R.R_Proj.transform_point(x+.009, y-.002, ccrs.Geodetic())
                        ax_n.text(trans[0], trans[1], lab, fontsize= 20)

        #####
        #Plot Radar Platforms (if desired and available)
        if isinstance(p, Radar):
            if p.type != 'MAINR':
                if print_long == True: print(p.name)
                #det if (any of) the radar is located within the area included in the plot domain at the time of the plot
                sites_subdf, valid_sites = p.grab_pform_subset(print_long, e_test, Data, bounding= t_R.Domain)
                
                #if these conditions are met then plot the radar marker(s)
                if valid_sites == True:
                    legend_elements.append(p.leg_entry)
                    if p.type == 'KA':
                        ## Plot the marker
                        ax_n.plot(p.lon, p.lat, transform=ccrs.PlateCarree(), marker=p.m_style, color=p.m_color, markersize=p.m_size,
                                  markeredgewidth=5, path_effects=[PathEffects.withStroke(linewidth=15, foreground='k')], zorder=10)
                        ## Plot RHI spokes
                        if config.rhi_ring == True: p.rhi_spokes_rings()
                        ## Optional textlabel on plot
                        if p.marker_label == True:
                            ax_n.text(p.lon+.009, p.lat-.002, p.name, transform=ccrs.PlateCarree(), path_effects=[PathEffects.withStroke(linewidth=4, foreground='xkcd:pale blue')])
                    if p.type == 'WSR':
                        ## Plot the marker
                        ax_n.plot(sites_subdf['lon'], sites_subdf['lat'], transform=ccrs.PlateCarree(), marker=p.m_style, color=p.m_color, markersize=p.m_size,
                                  markeredgewidth=5, path_effects=[PathEffects.withStroke(linewidth=12, foreground='k')], zorder=8)
                        # include labels for the sites on the plots
                        if p.marker_label == True:
                            for x, y, lab in zip(sites_subdf['lon'], sites_subdf['lat'], sites_subdf['R_Name']):
                                trans = t_R.R_Proj.transform_point(x, y, ccrs.Geodetic())
                                ax_n.text(trans[0]+.03, trans[1]-.01, lab, fontsize= 27)
                    if p.type == 'NOXP': print('Code not written yet')

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
    
    ## DEAL WITH COLORBARS
    # Attach colorbar to each subplot
    divider = make_axes_locatable(plt.gca())
    c_ax = divider.append_axes("right", "5%", pad="2%", axes_class=plt.Axes)
    sm = plt.cm.ScalarMappable(cmap=c_scale, norm=matplotlib.colors.Normalize(vmin=vminb, vmax=vmaxb))
    sm._A = []
    cb = plt.colorbar(sm, cax=c_ax, label=c_label)
    
    ## SET UP LEGENDS
    if leg == True: #this means you are currently making the left subplot
        #add legend for platform markers
        l = ax_n.legend(handles=legend_elements, loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(0,.5), 
                        handlelength=.1, borderpad=.5, title="Platforms", shadow=True, fancybox=True, ncol=1, edgecolor='black')
        l.get_title().set_fontweight('bold')
        l.get_title().set_fontsize(25)
    if leg == False:  #this means you are currently making the right subplot
        ## Plot platform colorbar
        #set up colorbar axis that will be as tall and 5% as wide as the 'parent' radar subplot
        cbar_ax = inset_axes(ax_n, width= '5%', height= '100%', loc='center left', bbox_transform=ax_n.transAxes, bbox_to_anchor=(-.25, 0,1,1)) 
        cbar = plt.colorbar(Data['Var'].CS3, cax=cbar_ax, orientation='vertical', label=Data['Var'].v_lab, ticks=MaxNLocator(integer=True))#,ticks=np.arange(Data['p_var'].global_min, Data['p_var'].global_max+1,2))

# * * * * * * *
def time_series(ts, Data, fig, sub_pos, print_long, e_test):
    ''' Plot a time series (sub)plot; could be of pvar info from various intstruments or of wind info 
    ----
    INPUTS:
    ts: ........fill in 
    Data: dict as described in ppiplot defn 
    fig & sub_pos: relating to the plot position within the larger figure (sub_pos is the gridspec info)
    print_long & e_test: bool str as described in the ppi defn 
    '''
    if print_long == True: print('~~~~~~~~~~~Made it into time_series~~~~~~~~~~~~~~~~~')

    ax_n= fig.add_subplot(sub_pos)

    TSleg_elements = [] #empty list to append the legend entries to for each subplot that is actually plotted
    ## MAKE THE TIMESERIES
    #### * * * * * * * * *
    if ts in ['Thetav', 'Thetae']:
        for p in Data.values():
            if isinstance(p, Torus_Insitu):
                if print_long == True: print('Plotting '+str(p.name)+' on time series')
                ## If the platform matches a type listed in TS_masking use masked data for the time series; else use the unmasked data
                if p.type in config.TS_masking[:]: plotting_data= p.ts_mask_df
                else: plotting_data= p.df[config.p_var].values

                ## Plot
                ax_n.plot(p.df['datetime'], plotting_data, linewidth=3, color=p.l_color) #assigning label= is what allows the legend to work
                TSleg_entry = Line2D([], [], label=p.leg_str, linewidth=12, color=p.l_color)
                TSleg_elements.append(TSleg_entry)

        ## Set up XY axes tick locations
        ax_n.set_ylim(bottom=Data['Var'].global_min, top=Data['Var'].global_max)
        ax_n.yaxis.set_major_locator(MultipleLocator(5)) # set up major tick marks (this is set up to go by 5's will want to change for diff vars)
        ax_n.yaxis.set_minor_locator(AutoMinorLocator(5)) # set up minor ticks (this have it increment by 1's will want to change for diff vars)
        ## Set up grid for plot
        ax_n.grid(which='major', axis='y', linestyle='--', linewidth=2)
        ax_n.grid(which='minor', axis='y', linestyle=':')
        ## Set up axes formats (lables/ ticks etc)
        ax_n.set_ylabel(Data['Var'].v_lab)

    # * * * 
    if ts == 'Wind':
        p = Data[config.Wind_Pform]
        if print_long == True: print('Plotting '+str(p.name)+' on time series')
        
        ax_n.plot(p.df['datetime'], p.df['spd'])
        ax_n.fill_between(p.df['datetime'], p.df['spd'], 0)
        TSleg_entry = Line2D([], [], label='Wind Spd', linewidth=12, color='tab:blue')
        TSleg_elements.append(TSleg_entry)

        ax2 = ax_n.twinx()
        ax2.plot(p.df['datetime'], p.df['dir'], '.k', linewidth=.05)
        TSleg_entry = Line2D([], [], marker='.', color='black', label='Wind Dir', markersize=26)
        TSleg_elements.append(TSleg_entry)

        ax_n.yaxis.set_major_locator(LinearLocator(numticks=5))
        ax_n.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        ax_n.axhline(0, color='k', linewidth=5, zorder=10)
        ax_n.set_ylabel('Wind Spd', multialignment='center')
        #set up minor ticks 
        #  ax2.yaxis.set_minor_locator(FixedLocator(np.arange(45, 405, 90)))
        #  ax2.grid(b=True, which='minor', axis='y', color='k', linestyle='--', linewidth=0.5)
        #  ax2.tick_params(which='minor', width=2, length=7, color='grey')
        #set up major tick marks 
        ax2.set_ylabel('Wind Dir ($^{\circ}$)', multialignment='center')
        ax2.set_ylim(0, 360)
        ax2.yaxis.set_major_locator(FixedLocator(np.arange(90, 360, 90)))
        ax2.grid(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=1)
        ax2.tick_params(which='major', width=2, length=14, color='black')

    '''

        ax2.yaxis.set_major_locator(FixedLocator(np.arange(90, 360, 90)))
        ax2.tick_params(which='major', width=2, length=14, color='black')
        ax2.grid(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=1)
        ax2.set_ylabel('Wind Dir ($^{\circ}$)', multialignment='center')
        ax2.set_ylim(0, 360)
    '''
    # * * * 
    ax_n.yaxis.set_label_coords(0, 0.1)
 
    # if desired this subsets the timerange that is displayed in the timeseries 
    if config.ts_extent != None:
        ax_n.set_xlim(Platform.Scan_time - timedelta(minutes=config.ts_extent), Platform.Scan_time + timedelta(minutes=config.ts_extent))

    #  ax_n.tick_params(axis='y', labelrotation=20)
    ax_n.margins(x=0, y=0)
    ax_n.xaxis.set_minor_locator(AutoMinorLocator(6)) # set up minor ticks (should be a multiple of ten intervals (ie 10,20,30... min spans)
    ax_n.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax_n.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M')) #strip the date from the datetime object
    ax_n.grid(which='minor', axis='x')
    ax_n.grid(which='major', axis='x', color='grey', linewidth=2.5)
    ax_n.set_xlabel('Time (UTC)')
    ax_n.tick_params(which='minor', width=2, length=7, color='grey')
    ax_n.tick_params(which='major', width=2, length=14, color='black')

    #if plotting more than one time series then only include the x axis label and ticks to the bottom timeseries
    num_of_TS= len(config.Time_Series)
    if num_of_TS != 1 and ax_n.rowNum != (num_of_TS-1):  
        plt.setp(ax_n.get_xticklabels(), visible=False) 

    ## If makeing the timeseries in conjunction with radar subplots set up vertical lines that indicate the time of
    #  radar scan and the timerange ploted (via filled colorline) on the radar plots
    if config.r_plotting == True:
        ax_n.axvline(Platform.Scan_time, color='r', linewidth=4, alpha=.5, zorder=10)
        ax_n.axvspan(Platform.Scan_time - timedelta(minutes=config.cline_extent), Platform.Scan_time + timedelta(minutes=config.cline_extent), facecolor='0.5', alpha=0.4, zorder=10)
    
    ## Include the legend
    leg = ax_n.legend(handles= TSleg_elements, loc='center left', shadow= True, fancybox= True, edgecolor='black')
    if ts == 'Wind':
        leg.set_title(config.Wind_Pform)
        leg.get_title().set_fontweight('bold')
        leg.get_title().set_fontsize(25)
    print(ax_n.lines)
    print('***')
    print(vars(ax_n))
    
    if print_long == True: print('~~~~~~~~~~~Made it through time_series~~~~~~~~~~~~~~')

##############################################################################################
#######################
## RADAR DEFINITIONS ##
#######################
def det_radar_fields(radar):
    #creating the mask for attenuation
    reflectivity = radar.fields['reflectivity']['data']
    spectrum_width = radar.fields['spectrum_width']['data']
    velocity = radar.fields['corrected_velocity']['data']
    total_power = radar.fields['total_power']['data']
    normal = radar.fields['normalized_coherent_power']['data']
    normal_mask = (normal.flatten() < 0.4)
    range_mask = np.zeros(np.shape(reflectivity))

    for i in range(0, len(range_mask[:,0])): range_mask[i,:] = radar.range['data'] > (radar.range['data'][-1]-1000.)

    range_mask = range_mask.astype(bool)
    total_mask = [any(t) for t in zip(range_mask.flatten(), normal_mask.flatten())]
    refl_mask = np.ma.MaskedArray(reflectivity, mask=normal_mask)
    sw_mask = np.ma.MaskedArray(spectrum_width, mask=normal_mask)
    vel_mask = np.ma.MaskedArray(velocity, mask=normal_mask)

    #create the dictionary for the masks
    refl_dict, sw_dict, vel_dict = {'data':refl_mask}, {'data':sw_mask}, {'data':vel_mask}
    radar.add_field('refl_fix', refl_dict)
    radar.add_field('sw_fix', sw_dict)
    radar.add_field('vel_fix', vel_dict)
# * * * 
def det_nearest_WSR(p_df):
    ''' locate the nearest WSR88D site to the specified insitu instruments
    '''
    #find the locations of all WSR88D sites (outputs dict in format {Site_ID:{lat:...,lon:...,elav:...], ...})
    all_WSR = pyart.io.nexrad_common.NEXRAD_LOCATIONS
    #  print(json.dumps(all_WSR_locs, sort_keys=True, indent=4))

    #set up empty dataframe with site Ids as column names
    d_from_all_r = pd.DataFrame(columns = all_WSR.keys())
    #fill in said dataframe with the distance from all 88D sites from each probe measurement
    for key in all_WSR:
        d_from_r = np.square(p_df['lat']-all_WSR[key]['lat']) + np.square(p_df['lon']-all_WSR[key]['lon'])
        d_from_all_r[key] = d_from_r
    #Determine which WS88D site is closest to the probe and add to the original probe dataframe
    p_df['Radar_ID'] = d_from_all_r.idxmin(axis = 1)
    
    #Determine the unique radar sites to be plotted 
    r_ofintrest = p_df.Radar_ID.unique()
    return r_ofintrest

# * * * 
def get_WSR_from_AWS(start, end, radar_id, download_directory):
    ''' Retrieve the NEXRAD files that fall within a timerange for a specified radar site from the AWS server
    ----------
    INPUTS start: datetime;  start of the desired timerange
           end: datetime; end of the desired timerange
           radar_id : string;  four letter radar designation
           download_directory: string; location for directory containing the downloaded radarfiles
    -------
    RETURN radar_list : Py-ART Radar Objects
    '''
    #Create this at the point of use otherwise it saves everything and eventually crashes 
    conn = nexradaws.NexradAwsInterface()

    #Det the radar scans that fall within the time randge for a given radar site
    scans = conn.get_avail_scans_in_range(start, end, radar_id)
    print("There are {} scans available between {} and {}\n".format(len(scans), start, end))
  
    #Dont download files that you alrady have ....
    path = config.temploc + config.day + '/radar/Nexrad/Nexrad_files'
    if not os.path.exists(path): Path(path).mkdir(parents=True)

    # missing_scans is a list of scans we don't have and need to download create_filepath returns tuple of 
    # (directory, directory+filename) [-1] returns the directory+filename
    missing_scans = list(filter(lambda x: not Path(x.create_filepath(path, False)[-1]).exists(), scans))

    # missing files is the list of filenames of files we need to download 
    missing_files = list(map(lambda x: x.create_filepath(path, False)[-1], missing_scans))
    print("missing "+ str(len(missing_files))+ " of "+ str(len(scans))+ " files\n"+ missing_files)

    results = conn.download(missing_scans, path, keep_aws_folders=False)

    print('{}\n{} downloads failed: {}\n'.format(results.success, results.failed_count, results.failed))
    #print("Results.iter_success : {}\n".format(reults.iter_success()))

    # missing_scans_after is a list of scans we don't have (download failed) create_filepath returns tuple of 
    # (directory, directory+filename) [-1] returns the directory+filename
    missing_files_after = list(filter(lambda x: not Path(x.create_filepath(path, False)[-1]).exists(), scans))

    if len(missing_files_after) > 0:
        print("ERROR: Some Radar Scans are Missing \n"+ str(missing_files_after))
        exit()

    radar_files = list(map(lambda x: x.create_filepath(path, False)[-1], scans))
    # Return list of files 
    return radar_files 

# * * * 
def read_from_nexrad_file(radar_file):
    radar = pyart.io.read_nexrad_archive(radar_file)
    return radar
def read_from_KA_file(radar_file):
    radar = pyart.io.read(radar_file)
    return radar
# Note: Cached version is cached on the file name, not the file contents.
# If file contents change you need to invalidate the cache or pass in the file contents directly to this function
#  function_cache_memory = Memory(config.g_cache_directory,verbose=1)
function_cache_memory = Memory(config.temploc, verbose=1)
cached_read_from_nexrad_file = function_cache_memory.cache( read_from_nexrad_file )
cached_read_from_KA_file = function_cache_memory.cache( read_from_KA_file )


# * * * 
def plot_radar_file(r_file, Data, subset_pnames, print_long, e_test, swp_id= None):
    print("open_pyart, scan file_name = {}\n".format(r_file))
    start_comptime = time.time()
    
    if config.Radar_Plot_Type == 'WSR_Plotting':
        #open file using pyart 
        try: radar = cached_read_from_nexrad_file(r_file)
        except: print("Failed to convert file: "+str(r_file))
    
    if config.Radar_Plot_Type == 'KA_Plotting':
        ## Read the radar file
        radar = cached_read_from_KA_file(r_file)
        ## If file contains a ppi scan proceed we are currently not interested in plotting the RHI scans
        if radar.scan_type == 'ppi':
            for swp_id in range(radar.nsweeps):
                tilt_ang = radar.get_elevation(swp_id) ## Det the actual tilt angle of a given sweep (returns an array)
                ## Check to see if the radarfile matches the elevation tilt we are interested in
                if np.around(tilt_ang[0], decimals=1) == config.p_tilt:
                    print("\nProducing Radar Plot:")
                    #  Assign radar fields and masking
                    det_radar_fields(radar)
    #  print(radar.info(level='compact'))

    ## Read in radar data and add to Data dict
    ##### + + + + + + + + + + + + + + + + + + +
    #  Establish info for the main plotting radar (aka scantime etc) & and locations for other radars (if deployed)
    Data, subset_pnames = Add_to_DATA('RADAR', Data, subset_pnames, print_long, MR_file=radar, swp=swp_id)
    if print_long == True: print(str(Data)+'\n')

    ## Proceed to plot the radar
    ##### + + + + + + + + + + + +
    ppiplot(Data, print_long, e_test, start_comptime)


##############################################################################################
#####################################################
# Read in Data that does not change for each image ##
#####################################################
print('\nRead in Torus Platforms')
#print(pform_names('ALL')) #list of ALL possible platforms
subset_pnames = [] #array of platforms that actually have data for the time range we are plotting
Data = {} #dictionary in which all the class objects will be stored (will contain platform data, locations, plotting variables etc)

## Read in the data for the TORUS Insitu platforms (if available)
Data, subset_pnames = Add_to_DATA('TInsitu', Data, subset_pnames, print_long)
## Establish info for the plotting variable (aka max min etc)
Data, subset_pnames = Add_to_DATA('PVar', Data, subset_pnames, print_long)

print('Read in Stationary Platforms Arrays\n')
## Read in the data for the Stationary Array platforms (if available)
Data, subset_pnames = Add_to_DATA('STN_I', Data, subset_pnames, print_long)


###################################
# Create Plot for each radarfile ##
###################################
## If Radar will be plotted
if config.r_plotting == True:
    print('\n Yes Plot Radar \n')
    # * * * 
    if config.Radar_Plot_Type == 'KA_Plotting':
        ## Get radar files
        radar_files = sorted(glob.glob(config.filesys+'TORUS_Data/'+config.day+'/radar/TTUKa/netcdf/*/dealiased_*'))
        ## Proceed to plot the radar
        ##### + + + + + + + + + + + +
        Parallel(n_jobs=config.nCPU, verbose=10)(delayed(plot_radar_file)(r_file, Data, subset_pnames, print_long, e_test) for r_file in radar_files)

    # * * * 
    if config.Radar_Plot_Type == 'WSR_Plotting': 
        #Det the unique radar sites to be plotted 
        unique_r_sites=det_nearest_WSR( Data[config.Centered_Pform].df)
        if print_long == True: print(unique_r_sites)

        #set up empty dataframe
        tranges_each_r = pd.DataFrame()
        for Rad_site in unique_r_sites:
            print(Rad_site)
            trange_r = Data[config.Centered_Pform].df.loc[Data[config.Centered_Pform].df.Radar_ID == Rad_site, ['datetime']].rename(columns={'datetime': Rad_site})
            trange_r_start, trange_r_end = trange_r.min(), trange_r.max()
            tranges_each_r = pd.concat([tranges_each_r, trange_r], axis=1)
            
            print("start "+str(trange_r_start[Rad_site])+ "\nend "+str(trange_r_end[Rad_site])+ "\n ***")
            radar_files = get_WSR_from_AWS(trange_r_start[Rad_site], trange_r_end[Rad_site], Rad_site, config.temploc)
            print('********\n Radar files to process:\n'+ str(radar_files))
            
            #Hard code the swp numbers that will be associated with a given tilt angle 
            if config.p_tilt == .5: swp_id=[0 , 1]
            elif config.p_tilt == 1: swp_id=[2 , 3]
            elif config.p_tilt == 1.5: swp_id=[4 , 5]
            else: print('The tilt angle {} is not hard coded yet for WSR'.format(config.p_tilt))

            #open the downloaded files as pyart objects
            print("Radar files to process:\n"+str(radar_files))

            ## Proceed to plot the radar
            ##### + + + + + + + + + + + +
            Parallel(n_jobs=config.nCPU, verbose=10)(delayed(plot_radar_file)(r_file, Data, subset_pnames, print_long, e_test, swp_id= swp_id) for r_file in radar_files)

        print(tranges_each_r)


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
print('\nIt took '+ str(time.time() - totalcompT_start)+ " to complete all of the plots\n")
print("ALL FINISHED")
