#import needed modules
######################
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patheffects as PathEffects
from matplotlib.ticker import (MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator)
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

## Imports form other files 
############################
import config #this is the file with the plotting controls to access any of the vars in that file use config.var

#rename a few commonly used vars so that the config.var does not have to be used repeatedly 
print_long, e_test = config.print_long, config.e_test

## Read in defns that I have stored in another file (for ease of use/consistancy accross multiple scripts)
from shared_defns import Add_to_DATA, pform_names, error_printing, Platform, Radar, Torus_Insitu, Stationary_Insitu 

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
    if print_long== True: print('~~~~~~~~~~~Made it into ppiplot~~~~~~~~~~~~~~~~~~~~~')

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
        fig = plt.figure(figsize=(20,10), facecolor='white')
        gs= GridSpec(nrows=2,ncols=1)
    else:
        ## Set up figure
        fig = plt.figure(figsize=(32,20),facecolor='white')
        ## Establish the gridspec layout
        gs= GridSpec(nrows=2, ncols=5,width_ratios=[.5,8,1,8,.25],height_ratios=[3,2],wspace=.25,hspace=.1)
        ## There are extra columns which are spacers to allow the formating to work
        # Uncomment to visualize the spacers
        #  ax_c=fig.add_subplot(gs[0,2])
        #  ax_l=fig.add_subplot(gs[0,0])
        #  ax_s=fig.add_subplot(gs[0,4])
    
    ## Plot each radar subplots #only set up for two radar subplots currently 
    col=-1
    for mom in config.r_mom[:]:
        # row and col pos of the subplot in the gridspec
          # row should always be topmost row (row= 0) 
          # first (left) subplot should be col=1 & the second (right) subplot should be col=3 due to the extra spacing "subplots" in the gridpec
        row, col= 0, col+2
        if print_long== True: print('Radar moment: ',mom,' Gridspec row: ',row,' GridSpec col: ',col)
        #if you are plotting the left subplot include the legend 
        if col == 1: leg= True 
        else: leg= False 
       
        ## Make the subplot
        #### * * * * * * * *
        post= radar_subplots(mom, Data, fig, gs[row,col], leg, print_long, e_test)
        
    ## Make the Times Series Subplots (if applicable)
    #### * * * * * * * * * * * * * * * * * ** * * * *
    if config.t_plotting == True: time_series(Data, fig, gs, print_long, e_test)

    ## Plot platform colorbar
    cbar_ax = plt.axes([.514, post.y0, .014, post.y1-post.y0])#left, bottom, width, height
    cbar=plt.colorbar(Data['Var'].CS3, cax=cbar_ax, orientation='vertical', label=Data['Var'].v_lab, ticks=MaxNLocator(integer=True))#,ticks=np.arange(Data['p_var'].global_min, Data['p_var'].global_max+1,2))

    ## Plot title
    plt.suptitle(Data['P_Radar'].name+' '+str(config.p_azimuth)+r'$^{\circ}$ PPI '+Data['P_Radar'].fancy_date_str, y=.92)

    ## Finish plot
    output_name=config.filesys+'TORUS_Data/'+config.day+'/mesonets/plots/'+Data['P_Radar'].name+'_'+config.p_var+'_'+Data['P_Radar'].time.strftime('%m%d_%H%M')+'.png'
    print(output_name)
    plt.savefig(output_name, bbox_inches='tight', pad_inches=.3)
    plt.close()

    ## Makes a ding noise
    print('\a')
    if print_long== True: print('~~~~~~~~~~~made it through ppiplot~~~~~~~~~~~~~~~~~~')
    print('Done Plotting \n ***************************************************************************************************')

# * * * * * *  *
def radar_subplots(mom, Data, fig, sub_pos, leg, print_long, e_test):
    ''' Plots each of the radar subplots and calls for the markers to be plotted for all the additional platforms
    ----
    INPUTS: mom: str, indicates which radar moment you would like to plot (ie Reflectivity, Velocity, Spectrum Width etc)
            fig, gs, & display: information about the plot
            ka_info, var, & pform: dictionarys, as described in the ppiplot comment
    '''
    if print_long== True: print('~~~~~~~~~~~made it into radar_subplots~~~~~~~~~~~~~~')
    
    ## SET UP VARS FOR EACH RADAR MOMENTS
    if mom == 'refl':
        p_title, field, c_label, c_scale, vminb, vmaxb= 'Reflectivity', 'refl_fix', 'Radar Reflectivity [dbz]', 'pyart_HomeyerRainbow', -30., 30.
    elif mom == 'vel':
        p_title, field, c_label, c_scale, vminb, vmaxb='Radial Velocity', 'vel_fix', 'Velocity [m/s]', 'pyart_balance', -40., 40.
    else: print("Hey what just happened!\n Check the Radar moments for spelling")

    #Assign the name of the Main Plotting Radar its own var
    #use Data[MRadar].var to find vars not stored in Data['P_Radar'] (like lat, lon etc) 
    MRadar=Data['P_Radar'].name

    ## SET UP SUBPLOTS
    ax_n=fig.add_subplot(sub_pos, projection=Data['P_Radar'].display.grid_projection)  
    ax_n.plot(Data[MRadar].lon, Data[MRadar].lat, transform=Data['P_Radar'].display.grid_projection)
    ax_n.text(.5, -.065, p_title, transform=ax_n.transAxes, horizontalalignment='center', fontsize=40) #the radar subplot titles
    Data['P_Radar'].display.plot_ppi_map(field, Data['P_Radar'].swp, title_flag=False, colorbar_flag=False, cmap=c_scale, ax=ax_n, vmin=vminb, vmax=vmaxb, embelish=False)

    ## PLOT PLATFORMS AS OVERLAYS(ie marker,colorline etc) ON RADAR
    #  iterate over each object contained in dict Data (returns the actual objects not their keys)
    legend_elements=[] #empty list to append the legend entries to for each platfrom that is actually plotted
    for p in Data.values():
        #print(vars(p))
        ###### 
        #Plot Inistu Torus Platforms (if desired and available)
        if isinstance(p, Torus_Insitu):
            # * * * 
            if config.NSSLm == True:
                if p.type == 'NSSL': 
                    if print_long== True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    p.plot_Tpform(Data, ax_n, print_long, e_test)
            # * * * 
            if config.NEBm == True:
                if p.type == 'UNL': 
                    if print_long== True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    p.plot_Tpform(Data, ax_n, print_long, e_test)
            # * * * 
            if config.UASm == True:
                if p.type == 'UAS': 
                    if print_long== True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    p.plot_Tpform(Data, ax_n, print_long, e_test, border_c='xkcd:very pale green', labelbias=(0,0.01))

        ###### 
        #Plot Stationary Inistu Platforms (if desired and available)
        if isinstance(p, Stationary_Insitu):
            # * * * 
            if config.MESONETSm == True:
                if p.type == 'WTM': 
                    if print_long== True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    ax_n.plot(Data['WTx_M'].df['lon'], Data['WTx_M'].df['lat'], data=Data['WTx_M'].df, transform=ccrs.PlateCarree(), marker=m_style, linestyle='None', markersize=15, color='k')
                    proj=display.grid_projection
                    #  WTM_sub=grab_platform_area_subset(WTM,ymin,ymax,xmin,xmax,p_var)
                    for x,y,lab in zip(Data['WTx_M'].df['lon'], Data['WTx_M'].df['lat'], Data['WTx_M'].df['Stn_ID']):
                        #  print(lab)
                        #  print(x)
                        #  print(y)
                        trans =proj.transform_point(x,y,ccrs.Geodetic())
                        ax_n.text(trans[0],trans[1],lab)
                    #  print(ax_n.transData.transform_point((31,-100)))
                    #  print(ax_n.get_xlim())
                    #  print(ax_n.get_ylim())
                    #  print(ax_n.transData.transform((30,100)))
                    #  print(ax_n.transData.transform_point((30,-100,ccrs.Geodetic())))
                    #  print(ax_n.projection.get_projection_class())
                    #      #ax_n.text(x,y,lab,transform=ax_n.transData)
                    #  text = [ax_n.text(*item) for item in zip(WTM['Lon'], WTM['Lat'], WTM['Stn_ID'])]
                    #ax_n.text(WTM['Lon'].all(),WTM['Lat'].all(), WTM['Stn_ID'].all(),transform=ccrs.PlateCarree()) 
                    #ax_n.plot(WTM['Lon'],WTM['Lat'], transform=ccrs.PlateCarree(), marker=r'$\AA$',markersize=10, color='k',zorder=9) 
            # * * * 
            if config.METARm == True:
                if p.type == 'METAR': 
                    if print_long== True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    print("To Be filled in ")
            # * * * 
            if config.ASOS_AWOSm == True:
                if p.type == 'ASOS/AWOS': 
                    if print_long== True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    print("To Be filled in ")

        ##### 
        #Plot Radar Platforms (if desired and available)
        if isinstance(p, Radar):
            # * * *
            if config.KAm == True: 
                if p.type == 'KA':
                    if print_long== True: print(p.name)
                    #if the radar is located within the area included in the plot
                    p_deploy = p.grab_pform_subset(print_long, e_test, Data, bounding= Platform.P_domain)
                    if p_deploy == True:
                        legend_elements.append(p.leg_entry)
                        ## Plot the marker
                        ax_n.plot(p.lon, p.lat, marker=p.m_style, transform=ccrs.PlateCarree(), color=p.m_color, markersize=18, 
                                  markeredgewidth=5, path_effects=[PathEffects.withStroke(linewidth=15,foreground='k')],zorder=10)
                        ## Plot RHI spokes 
                        if config.rhi_ring == True: p.rhi_spokes_rings() 
                        ## Optional textlabel on plot
                        #    d2=ax_n.text(p.lon+.006, p.lat-.011, p.name, transform=ccrs.PlateCarree(), path_effects=[PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')])
            # * * * 
            if config.WSRm == True: 
                if p.type == 'WSR':
                    if print_long== True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    ## Plot the marker
                    ax_n.plot(p.df['lon'], p.df['lat'], data=p.df, transform=ccrs.PlateCarree(), marker=m_style, linestyle='None', markersize=15, color='k')
                    #  ax_n.plot(p.lon, p.lat, marker=p.m_style, transform=ccrs.PlateCarree(), color=p.m_color, markersize=18,
                              #  markeredgewidth=5, path_effects=[PathEffects.withStroke(linewidth=15,foreground='k')],zorder=10)
            # * * * 
            if config.NOXPm == True: 
                if p.type == 'NOXP':
                    if print_long== True: print(p.name)
                    legend_elements.append(p.leg_entry)
                    print('Code not written yet')

    ## DEAL WITH COLORBARS
    # Attach colorbar to each subplot
    divider = make_axes_locatable(plt.gca())
    c_ax = divider.append_axes("right", "5%", pad="2%", axes_class=plt.Axes)
    sm = plt.cm.ScalarMappable(cmap=c_scale, norm=matplotlib.colors.Normalize(vmin=vminb, vmax=vmaxb))
    sm._A = []
    cb = plt.colorbar(sm, cax=c_ax, label=c_label)
    if leg== False:  #this means you are currently making the right subplot
        post= ax_n.get_position() ## Get position information to pass along to the remaining colorbar

    ## SET UP LEGENDS
    if leg == True: #this means you are currently making the left subplot
        #add legend for platform markers
        l=ax_n.legend(handles=legend_elements,loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(0,.5), handlelength=.1,title="Platforms",shadow=True,fancybox=True,ncol=1,edgecolor='black')
        l.get_title().set_fontweight('bold')
        post= np.nan

    ## PLOT BACKGROUND FEATURES
    if config.country_roads == True:
        ox.config(log_file=True, log_console=True, use_cache=True) #the config in this line has nothing to do with config.py
        G = ox.graph_from_bbox(Platform.P_domain.ymax, Platform.P_domain.ymin, Platform.P_domain.xmax, Platform.P_domain.xmin)
        ox.save_load.save_graph_shapefile(G, filename='tmp'+str(0), folder=config.filesys+'radarproc/roads/', encoding='utf-8')
        fname = config.filesys+'radarproc/roads/tmp'+str(0)+'/edges/edges.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree(), edgecolor='gray', linewidth=0.5)
        ax_n.add_feature(shape_feature, facecolor='none')
        shutil.rmtree(config.filesys+'radarproc/roads/tmp'+str(0)+'/')
    if config.hwys == True:
        fname = config.filesys+'radarproc/roads/GPhighways.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree(), edgecolor='grey')#edgecolor='black')
        ax_n.add_feature(shape_feature, facecolor='none')
    if config.county_lines == True:
        fname = config.filesys+'radarproc/roads/cb_2017_us_county_5m.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree(), edgecolor='gray')
        ax_n.add_feature(shape_feature, facecolor='none', linewidth=1.5, linestyle="--")
    if config.state_lines == True:
        states_provinces = cartopy.feature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
        ax_n.add_feature(states_provinces, edgecolor='black', linewidth=2)
    if print_long== True: print('~~~~~~~~~~~Made it through radar_subplots~~~~~~~~~~~')
    return post

# * * * * * * *
def plot_winds(self, ws, wd, wsmax, plot_range=None):
        """
        Required input:
            ws: Wind speeds (knots)
            wd: Wind direction (degrees)
            wsmax: Wind gust (knots)
        Optional Input:
            plot_range: Data range for making figure (list of (min,max,step))
        """ 
        self.start = dates[0]
        self.end = dates[-1]
        self.axis_num = 0
        self.dates = mpl.dates.date2num(dates)
        self.time = time.strftime('%Y-%m-%d %H:%M UTC')

                # set the starttime and endtime for plotting, 24 hour range
        endtime = dt.datetime(2016, 3, 31, 22, 0, 0, 0)
        starttime = endtime - dt.timedelta(hours=24)


        # Parse dates from .csv file, knowing their format as a string and convert to datetime
        def parse_date(date):
            return dt.datetime.strptime(date.decode('ascii'), '%Y-%m-%d %H:%M:%S')


        testdata = np.genfromtxt(get_test_data('timeseries.csv', False), names=True, dtype=None,
                                 usecols=list(range(1, 8)),
                                 converters={'DATE': parse_date}, delimiter=',')

        # Temporary variables for ease
        ws = testdata['WS']
        wsmax = testdata['WSMAX']
        wd = testdata['WD']
        date = testdata['DATE']

        data = {'wind_speed': (np.array(ws) * units('m/s')).to(units('knots')),
                'wind_speed_max': (np.array(wsmax) * units('m/s')).to(units('knots')),
                'wind_direction': np.array(wd) * units('degrees')}

        meteogram = Meteogram(fig, data['times'], probe_id)
        meteogram.plot_winds(data['wind_speed'], data['wind_direction'], data['wind_speed_max'])

        ############################################

        # PLOT WIND SPEED AND WIND DIRECTION
        ln1 = self.ax_n.plot(self.dates, ws, label='Wind Speed')
        self.ax_n.fill_between(self.dates, ws, 0)
        self.ax_n.set_xlim(self.start, self.end)
        if not plot_range: plot_range = [0, 20, 1]
        self.ax_n.set_ylabel('Wind Speed (knots)', multialignment='center')
        self.ax_n.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        self.ax_n.grid(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=0.5)
        ln2 = self.ax_n.plot(self.dates, wsmax, '.r', label='3-sec Wind Speed Max')

        ax7 = self.ax_n.twinx()
        ln3 = ax7.plot(self.dates, wd, '.k', linewidth=0.5, label='Wind Direction')
        ax7.set_ylabel('Wind\nDirection\n(degrees)', multialignment='center')
        ax7.set_ylim(0, 360)
        ax7.set_yticks(np.arange(45, 405, 90), ['NE', 'SE', 'SW', 'NW'])
        lines = ln1 + ln2 + ln3
        labs = [line.get_label() for line in lines]
        ax7.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H UTC'))
        ax7.legend(lines, labs, loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3, prop={'size': 12})

# * * * * * * *
def time_series(Data, fig, gs, print_long, e_test):
    ''' Plot the time series of p_var from the various instrument platforms
    ----
    INPUTS; fig & gs: relating to the plot layout
            ka_info, var, & pform: dictionarys, as described in the ppiplot comment
            r_plotting: True/False, If False will only plot timeseries (no radar).. In theroy have not actually done this yet
    '''
    if print_long== True: print('~~~~~~~~~~~Made it into time_series~~~~~~~~~~~~~~~~~')

    if config.r_plotting == True: ax_n=fig.add_subplot(gs[1,:])
    else: fig, ax_n= plt.subplots() #will only plot the time series (not the radar subplots)

    ## MAKE THE TIMESERIES
    #### * * * * * * * * *
    for p in Data.values():
        if isinstance(p, Torus_Insitu):
            if print_long== True: print('Plotting ',p.name,' on time series')
            ## If the platform matches a type listed in TS_masking use masked data for the time series; else use the unmasked data
            if p.type in config.TS_masking[:]: plotting_data=p.ts_mask_df
            else: plotting_data= p.df[config.p_var].values

            wind=True
            ## Plot
            if wind == False:
                ax_n.plot(p.df['datetime'], plotting_data, linewidth=3, color=p.l_color, label=p.leg_str) #assigning label= is what allows the legend to work

            if wind == True:
                p=Data['Prb1']
                l1= ax_n.plot(p.df['datetime'], p.df['spd'], label='Wind Speed')
                ax_n.fill_between(p.df['datetime'], p.df['spd'], 0)
                ax_n.set_xlim(config.tstart, config.tend)
                ax_n.set_ylabel('Wind Speed (knots)', multialignment='center')

                l2= ax_n.plot(p.df['datetime'], p.df['spd'], p.df['wmax'],'.r', label='30-sec Wind Speed Max')
                
                ax7 = ax_n.twinx()
                l3= ax7.plot(p.df['datetime'], p.df['dir'], '.k', linewidth=.5, label='Wind Direction')
                ax7.set_ylabel('Wind\nDirection\n(degrees)', multialignment='center')
                ax7.set_ylim(0, 360)
                ax7.set_yticks(np.arange(45, 405, 90), ['NE', 'SE', 'SW', 'NW'])
                lines = ln1 + ln2 + ln3
                labs = [line.get_label() for line in lines]
                ax7.legend(lines, labs, loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3, prop={'size': 12})
    
    ## Set up XY axes tick locations
    ax_n.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M')) #strip the date from the datetime object
    ax_n.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax_n.xaxis.set_minor_locator(AutoMinorLocator(6)) # set up minor ticks (should be a multiple of ten intervals (ie 10,20,30... min spans)
    ax_n.yaxis.set_major_locator(MultipleLocator(5)) # set up major tick marks (this is set up to go by 5's will want to change for diff vars)
    ax_n.yaxis.set_minor_locator(AutoMinorLocator(5)) # set up minor ticks (this have it increment by 1's will want to change for diff vars)

    ## Set up axes formats (lables/ ticks etc)
    ax_n.set_ylabel(Data['Var'].v_lab)
    ax_n.set_xlabel('Time (UTC)')
    ax_n.tick_params(which='major', width=2,length=20, color='black')
    ax_n.tick_params(which='minor', width=2,length=10, color='grey')

    ## Set up grid for plot
    ax_n.grid(True)
    ax_n.grid(which='major', axis='x', color='grey', linewidth=2.5)
    ax_n.grid(which='minor', axis='x')
    ax_n.grid(which='major', axis='y', linestyle='--', linewidth=2)
    ax_n.grid(which='minor', axis='y', linestyle=':')

    ## If makeing the timeseries in conjunction with radar subplots set up vertical lines that indicate the time of 
    #  radar scan and the timerange ploted (via filled colorline) on the radar plots
    if config.r_plotting == True:
        ax_n.axvline(Data['P_Radar'].time, color='r', linewidth=4, alpha=.5, zorder=1)
        ax_n.axvspan(Data['P_Radar'].time-timedelta(minutes=config.cline_extent), Data['P_Radar'].time+timedelta(minutes=config.cline_extent), facecolor='0.5', alpha=0.4)

    ## Include the legend
    leg=ax_n.legend()
    for line in leg.get_lines(): line.set_linewidth(12)

    if print_long== True: print('~~~~~~~~~~~Made it through time_series~~~~~~~~~~~~~~')

##############################################################################################
#######################
## RADAR DEFINITIONS ##
#######################
def det_radar_feilds(radar):
    radar_name = radar.metadata['instrument_name']

    #creating the mask for attenuation
    reflectivity = radar.fields['reflectivity']['data']
    spectrum_width = radar.fields['spectrum_width']['data']
    velocity = radar.fields['corrected_velocity']['data']
    total_power = radar.fields['total_power']['data']
    normal = radar.fields['normalized_coherent_power']['data']
    normal_mask = (normal.flatten() < 0.4)
    range_mask = np.zeros(np.shape(reflectivity))

    for i in range(0,len(range_mask[:,0])): range_mask[i,:] = radar.range['data']>(radar.range['data'][-1]-1000.)

    range_mask = range_mask.astype(bool)
    total_mask = [any(t) for t in zip(range_mask.flatten(), normal_mask.flatten())]
    refl_mask = np.ma.MaskedArray(reflectivity, mask=normal_mask)
    sw_mask = np.ma.MaskedArray(spectrum_width, mask=normal_mask)
    vel_mask = np.ma.MaskedArray(velocity, mask=normal_mask)

    #create the dictionary for the masks
    refl_dict = {'data':refl_mask}
    sw_dict = {'data':sw_mask}
    vel_dict = {'data':vel_mask}
    radar.add_field('refl_fix',refl_dict)
    radar.add_field('sw_fix',sw_dict)
    radar.add_field('vel_fix',vel_dict)

##############################################################################################
#####################################################
# Read in Data that does not change for each image ##
#####################################################

#print(pform_names('ALL')) #list of ALL possible platforms 
subset_pnames=[] #array of platforms that actually have data for the time range we are plotting 
Data={} #dictionary in which all the class objects will be stored (will contain platform data, locations, plotting variables etc)

## Read in the data for the TORUS Insitu platforms (if available)
Data, subset_pnames = Add_to_DATA('TInsitu',Data,subset_pnames,print_long)
## Establish info for the plotting variable (aka max min etc) 
Data, subset_pnames = Add_to_DATA('PVar',Data,subset_pnames,print_long)


###################################
# Create Plot for each radarfile ##
###################################
## If Radar will be plotted 
if config.r_plotting == True:
    print('\n Yes Print Radar \n')
    ## Get radar files
    fil =sorted(glob.glob(config.filesys+'TORUS_Data/'+config.day+'/radar/TTUKa/netcdf/*/dealiased_*'))
    #  ka1_files =sorted(glob.glob(config.filesys+'TORUS_Data/'+config.day+'/radar/TTUKa/netcdf/ka1/dealiased_*'))
    #  ka2_files =sorted(glob.glob(config.filesys+'TORUS_Data/'+config.day+'/radar/TTUKa/netcdf/ka2/dealiased_*'))
    #  fil= ka1_files + ka2_files

    ## Read in and check each radarfile in fil 
    #### + + + + + + + + + + + + + + + + + + +
    for thefile in fil[:]:
        print(str(thefile))
        radar = pyart.io.read(thefile) ## Read the radar file

        ## We are currently not interested in plotting the RHI scans 
        if radar.scan_type == 'rhi': print("This scan is an RHI \n **********************************")
        ## If file contains a ppi scan proceed 
        elif radar.scan_type == 'ppi':
            n_swps = radar.nsweeps
            for swp_id in range(n_swps):
                plotter = pyart.graph.RadarDisplay(radar)
                azimuth = radar.fixed_angle['data'][swp_id] #det the actual azimuth of a given sweep 

                ## Check to see if the radarfile matches the azimuth we are interested in 
                if np.around(azimuth, decimals=1) == config.p_azimuth: 
                    print("Producing Radar Plot")
                    det_radar_feilds(radar)  #Assign radar feilds and masking

                    ## Det Main plotting radar (aka which radar is associated with rfile) 
                    ##### + + +  + + + + + + + + + + + + + + + + + + + + + + + + + + + +
                    if thefile.split("/")[-1][10:13] == 'Ka1': P_Rname = 'Ka1'
                    elif thefile.split("/")[-1][10:13] == 'Ka2': P_Rname = 'Ka2'
                    else: print('What radar are we using?')

                    ## Read in radar data and add to Data dict
                    ##### + + + + + + + + + + + + + + + + + + +
                    #  Establish info for the main plotting radar (aka scantime, azimuth, name, loc, etc)
                    Data, subset_pnames = Add_to_DATA('Plotting_Radar',Data,subset_pnames,print_long, rfile=radar, rname=P_Rname, swp=swp_id)
                    #  Define the scantime for the radarfile 
                    Platform.Scan_time = Data['P_Radar'].time 
                    #  Define the plotting domain (centered on the plotting radar in this case)
                    Platform.P_domain = Data[Data['P_Radar'].name].getLocation(offsetkm=config.offsetkm)
                    #  Define pyart display object for plotting radarfile  
                    Radar.display = pyart.graph.RadarMapDisplay(Data['P_Radar'].rfile) 
                    #  set the limits in the display object (this will enforce the desired plotting domain)
                    Data['P_Radar'].display.set_limits(xlim=(Platform.P_domain.xmin, Platform.P_domain.xmax), ylim=(Platform.P_domain.ymin, Platform.P_domain.ymax))    
                    #  Establish locations for other radars (if deployed)
                    Data, subset_pnames = Add_to_DATA('RADAR', Data, subset_pnames, print_long)
       
                    ## Proceed to plot the radar
                    ##### + + + + + + + + + + + +
                    ppiplot(Data, print_long, e_test)

                ## If the scans azimuth does not match the one we are plotting
                else: print("Scan does not match the azimuth currently being plotted \n Currently plotting: azimuth= " 
                             + str(config.p_azimuth) +"\n This scan: azimuth= "+ str(azimuth)+'\n**********************************')


################################
# Create Timeseries only plot ##
################################
#Only plot timeseries (this code isn't fully fleshed out but in theroy this code is built in such a way to allow for this)
if config.r_plotting == False and config.t_plotting== True:
    print("Plot Timeseries only \n"+ str(filesys+'TORUS_Data/'+day+'/mesonets/NSSL/*.nc'))
    time_series(Data,r_plotting=False)
    fig.savefig('test2.png')
    plt.close()

###########################
print("ALL FINISHED")
