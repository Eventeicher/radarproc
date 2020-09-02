import numpy as np
from datetime import datetime
from datetime import timedelta
import os
import tempfile
from boto.s3.connection import S3Connection
import matplotlib.patheffects as mpatheffects
import matplotlib.pyplot as plt
from netCDF4 import num2date
import pyart
import cmocean
from matplotlib.ticker import MaxNLocator
import pandas as pd
import glob
import matplotlib
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from metpy.plots import ctables
import metpy
import matplotlib.patheffects as patheffects
import cartopy
import cartopy.crs as ccrs
import datetime as dt
import shapely.geometry as sgeom
from pyproj import Geod
import pandas as pd
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
import sys, traceback
from read_platforms import read_nsslmm, read_unlmm, maskdata
import json
import nexradaws
import pytz
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cftime

## Imports form other files
############################
import config #this is the file with the plotting controls to access any of the vars in that file use config.var

#rename a few commonly used vars so that the config.var does not have to be used repeatedly
print_long, e_test = config.print_long, config.e_test

## Read in defns that I have stored in another file (for ease of use/consistancy accross multiple scripts)
from shared_defns import Add_to_DATA, pform_names, error_printing, Platform, Radar, R_Plt, Torus_Insitu, Stationary_Insitu, get_WSR_from_AWS


from pathlib import Path
from os.path import expanduser

from joblib import Memory

from joblib import Parallel, delayed

# To run with a) warnings and b) stack trace on abort
# python3 -Walways  -q -X faulthandler plot_nexrad_insitu.py

##########################################################################
## VARIABLES
############
#day='20190524' #'YYYYMMDD'
day='20190517' #'YYYYMMDD'    # CV - These are the only files I have...

radar_plotting= True #would you like to plot radar images? (set to False to only plot timeseries)
r_only= False #Set to True for only radar as output, Set to False for Radar + Timeseries
p_var = "Thetav" #which var to plot (current options; Thetae, Thetav)
probe_of_interest='Prb1'

function_cache_memory = Memory(config.g_cache_directory,verbose=1)

#Troubleshooting y/n
####################
#Set to True to print statements when entering/leaving definitions
print_long = False # True/False variable

#Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting
  ##there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = False# True/False variable


##########################################################################
#################################
## KA & NOXP RADAR DEFENITIONS ##
#################################
def det_radar_deps(currentscantime, r_testing, e_test):
    ''' Determine if a given (ka) radar is deployed and if so assign the correct values to it.
        r_testing is the name of the radar you are testing to see if deployed
    '''
    #test if radar deployed
    try:
        kadep=pd.read_csv(config.g_root + 'TORUS_Data/'+day+'/radar/TTUKa/csv/'+day+'_deployments_'+r_testing+'.csv')

        #get radar positioning data
        k_lat, k_lon, head_ing, rhi_b, rhi_e = getRadarPositionData(currentscantime, kadep)

    except:
        error_printing(e_test)
        #If one of the radars did not deploy at all on a given day this will set its values to np.nan
        kadep, k_lat, k_lon, head_ing, rhi_b, rhi_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    return kadep, k_lat, k_lon, head_ing, rhi_b, rhi_e

# * * * * * * *
def getRadarPositionData(c_time, r_dep):
    ''' If a Ka Radar had deployed this funcion returns it's location (lat/lon) and
    the first and last RHI scan angle of a sweep if applicable.
    '''
    define_check=0
    for t in range(r_dep.time_begin.count()):
        beginscan=datetime.strptime(r_dep.time_begin[t], "%m/%d/%Y %H:%M")
        endscan=datetime.strptime(r_dep.time_end[t], "%m/%d/%Y %H:%M")

        if c_time >= beginscan and c_time <= endscan:
            try:
                klat, klon, head = r_dep.lat[t], r_dep.lon[t], r_dep.heading[t]
            except:
                klat, klon, head = np.nan, np.nan, np.nan
                error_printing(e_test)

            try:
                rhib, rhie = r_dep.rhib[t], r_dep.rhie[t]
            except:
                rhib, rhie = np.nan, np.nan
                error_printing(e_test)
        else:
            rhib, rhie, klat, klon, head = np.nan, np.nan, np.nan, np.nan, np.nan

        # preseve the actual values so they are not overwritten by other deployments
        if str(klat) != 'nan':
            KLAT, KLON, HEAD, RHIB, RHIE = klat, klon, head, rhib, rhie
            define_check = 1
        else: pass

    #In the instance that the second radar isnt deployed at time of radarscan the KLAT, KLON,.... will still be defined
      ## If the other radar is deployed it will not overwrite the assigned values
    if define_check == 0:
        KLAT, KLON, HEAD, RHIB, RHIE = np.nan, np.nan, np.nan, np.nan, np.nan
    else: pass

    return KLAT, KLON, HEAD, RHIB, RHIE



#################################
## RADAR RETRIEVAL DEFENITIONS ##
#################################
def det_nearest_WSR(p_df):
    '''
    Function to locate the nearest WSR88D site to the insitu instruments
    '''
    #find the locations of all WSR88D sites (outputs dictionary in format {Site_ID:{lat:...,lon:...,elav:...], ...}
    locs = pyart.io.nexrad_common.NEXRAD_LOCATIONS
    #print(json.dumps(locs, sort_keys=True, indent=4))

    #set up empty dataframe with site Ids as column names
    d_from_all_r = pd.DataFrame(columns=locs.keys())

    #fill in dataframe with the distance from all 88D sites from each probe measurement
    for key in locs:
        d_from_r=np.square(p_df['lat']-locs[key]['lat'])+np.square(p_df['lon']-locs[key]['lon'])
        d_from_all_r[key]=d_from_r

    #Determine which WS88D site is closest to the probe and add to the original probe dataframe
    p_df['Radar_ID']=d_from_all_r.idxmin(axis=1)
    #print(p_df)

    return p_df


###########################
## PLOTTING DEFINTIONS  ###
###########################
def getLocation(file, currentscantime, e_test, offsetkm, given_bearing=False):
    '''
    This definition has two functions:
        1) If no bearing is specified it will return a dict containing the max/min lat/lons
            to form a square surrounding the point indicated by lat1,lon1 by x km.
        2) If a bearing is given then the defintion will return one set of lat/lon values
             (end_lat and end_lon) which is the location x km away from the point lat1, lon1
             if an observer moved in a straight line the direction the bearing indicated.
    ----
    INPUTS
    file, currentscantime and var: used to find current lat and lon
    current_lon & current_lat: the starting position
    offsetkm: distance traveled from the starting point to the new locations
    given_bearing: True/False, are you given a specified direction of "travel"
    '''

    #determine the position of the probe at the time of radar scan
    df_sub, p_deploy = cached_grab_platform_subset(file, print_long, e_test, scan_time=currentscantime)

    current_lon=df_sub['lon'][int(len(df_sub['lon'])/2)]
    current_lat=df_sub['lat'][int(len(df_sub['lat'])/2)]

    #    error_printing(e_test)


    lat1 = current_lat * np.pi / 180.0
    lon1 = current_lon * np.pi / 180.0

    R = 6378.1 #earth radius
    #R = ~ 3959 MilesR = 3959

    if given_bearing == False:
        for brng in [0,90,180,270]:
            bearing = (brng / 90.)* np.pi / 2.

            new_lat = np.arcsin(np.sin(lat1) * np.cos(offsetkm/R) + np.cos(lat1) * np.sin(offsetkm/R) * np.cos(bearing))
            new_lon = lon1 + np.arctan2(np.sin(bearing)*np.sin(offsetkm/R)*np.cos(lat1),np.cos(offsetkm/R)-np.sin(lat1)*np.sin(new_lat))
            new_lon = 180.0 * new_lon / np.pi
            new_lat = 180.0 * new_lat / np.pi

            if brng == 0: max_lat=new_lat
            elif brng == 90: max_lon= new_lon
            elif brng == 180: min_lat= new_lat
            elif brng == 270: min_lon=new_lon

        bndry={'ymin': min_lat, 'ymax': max_lat, 'xmin': min_lon, 'xmax': max_lon}
        #return both individual points and dict for ease of use [same info but sometimes one format is easier than the other]
        return min_lat, max_lat, min_lon, max_lon, bndry

    else:
        #if a bearing is provided
        bearing = (given_bearing/ 90.)* np.pi / 2.

        new_lat = np.arcsin(np.sin(lat1) * np.cos(offsetkm/R) + np.cos(lat1) * np.sin(offsetkm/R) * np.cos(bearing))
        new_lon = lon1 + np.arctan2(np.sin(bearing)*np.sin(offsetkm/R)*np.cos(lat1),np.cos(offsetkm/R)-np.sin(lat1)*np.sin(new_lat))
        end_lon = 180.0 * new_lon / np.pi
        end_lat = 180.0 * new_lat / np.pi

        return end_lat, end_lon

#* * * * *
# This was removed but still seems to be needed???
# So I brought it back to get it to complie.
#

def legend_maker(p_name,m_style,m_color,leg_str,legend_elements):
    if print_long== True: print('made it into legend_maker')

    if (p_name in ['Ka1','Ka2']): #the legend entries for the KA radars
        legend_entry=Line2D([], [], marker=m_style, markeredgecolor='black',markeredgewidth=3,label=leg_str,markerfacecolor=m_color, markersize=26)
    else:
        legend_entry=Line2D([], [], marker=m_style, markeredgecolor=m_color,markeredgewidth=3,label=leg_str, markersize=26,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')])

    legend_elements=np.append(legend_elements,legend_entry)

    if print_long== True: print('made it through legend maker')
    return legend_elements


# This was removed but still seems to be needed???
# So I brought it back to get it to complie.

def platform_attr(file_type,legend_elements,radar_m=False,r_s=False):
    if print_long== True: print('Made it into platform_attr')

    print("file_type = ", file_type)

    m_style, m_color, l_color, leg_str= '1','red','red','FixMe! There was not a case for this file type'

    #assign the atributes for non-radar markers
    if radar_m == False:
        if file_type == "Prb1":
            m_style, m_color, l_color, leg_str= '1','xkcd:lightblue','steelblue','Prb1'
        elif file_type == "Prb2":
            m_style, m_color, l_color, leg_str= '1','xkcd:watermelon','xkcd:dusty red','Prb2'
        elif file_type == "FFld":
            m_style, m_color, l_color, leg_str= '1','xkcd:bubblegum pink','xkcd:pig pink','FFld'
        elif file_type == "LIDR":
            m_style, m_color, l_color, leg_str= '1','xkcd:pastel purple','xkcd:light plum','LIDR'
        elif file_type == "WinS":
            m_style, m_color, l_color, leg_str= '1','xkcd:peach','xkcd:dark peach','WinS'
        elif file_type == "CoMeT1":
            m_style, m_color, l_color, leg_str= '1','brown','brown','CoMeT1'
        elif file_type == "CoMeT2":
            m_style, m_color, l_color, leg_str= '1','yellow','yellow','CoMeT2'
        elif file_type == "CoMeT3":
            m_style, m_color, l_color, leg_str= '1','black','black','CoMeT3'

    #assign the atributes for the radar markers
    elif radar_m == True:
        if file_type == 'Ka2':
            m_style, m_color, l_color, leg_str= '8','xkcd:crimson','xkcd:crimson','Ka2'
        elif file_type == 'Ka1':
            m_style, m_color, l_color, leg_str= '8','mediumseagreen','mediumseagreen','Ka1'

    if r_s == True: #only add to the legend if the platform_attr def was called while making the radar subplots
        legend_elements=legend_maker(file_type, m_style,m_color,leg_str,legend_elements)

    if print_long== True: print('Made it through platform_attr')

    p_attr = {'m_style': m_style, 'm_color': m_color, 'l_color': l_color, 'leg_str': leg_str}

    # return m_style, m_color, l_color, leg_str, legend_elements

    return p_attr, legend_elements

# * * * * * *
def plot_colourline(x,y,c,cmap,ax,datacrs,color,amin=None,amax=None):
    ''' This defintion plots line extending out from the platform marker +/- in time
            The colorfill indicates values of the specifiec p_var (ie Thetae etc)
        ---
        INPUTS
        x & y: the loction of the probe at a given time
        c: the value of p_var at a given time (will correspond eith the colorfill)
        cmap: the colormap for the colorfill,
        ax & datacrs: the subplot to plot the line on, the projection of the data
        color: color of the border
        amin & amax: max and min value of p_var. prevents outliers from sckewing the color scale.
    '''
    if amin: pass
    else: amin=np.nanmin(c)

    if amax: pass
    else: amax=np.nanmax(c)

    c = cmap((c-amin)/(amax-amin))
    ax = plt.gca()
    for i in np.arange(len(x)-1): #This is the border of the colorline
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=color, linewidth=10.5, transform=datacrs)
    for i in np.arange(len(x)-1): #This is the colorramp colorline
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=c[i], linewidth=7.5, transform=datacrs)

    return

# * * * * * *
def grab_platform_subset(p_df, print_long, e_test, bounding=None, scan_time=None, time_offset= 300):
    ''' This def will take a given dataset and subset it either spatially or temporially
            1) If scan_time is given the data will be subset temporally
                    Grabs the observed thermo data +/- x seconds around radar scan_time.
                    The length of the subset is det by time_offset which is defaults to +/- 300 secs.
                    This means that unless a different time_offeset is specified the subset will span 10 min.
                    NOTE: This does not automatically control the gray vertical box on the timeseries
            2) If bounding is given the data will be subset spacially
                    will return the data points that fall within the box defined by ymin,ymax,xmin,xmax
            3) If both scan_time and a bounding region are provided the dataset should be subset both
                    temporally and spatially. (This has not been tested yet so double check if using)
    -----
    Returns: The subset dataset (df_sub) and a True/False statement regarding any data in the original dataset
                matched the subsetting criteria (p_deploy)
    '''

    #Temporal subset
    if scan_time != None:
        #aaa = p_df.loc[(p_df['datetime']>=scan_time-dt.timedelta(seconds=time_offset))]
        #df_sub = aaa.loc[(aaa['datetime']<=scan_time+dt.timedelta(seconds=time_offset))]
        #if print_long == True: print('Dataset has been temporally subset')


        # Convert scan_time from netcdf cftime to python standard datetime
        scan_time = pd.to_datetime(scan_time.strftime())

        td = dt.timedelta(seconds=time_offset)
        start = scan_time - td;
        end = scan_time + td;

        # Extract scans that fall in time range
        p_df['datetime'] = pd.to_datetime(p_df['datetime']) # Normalize time to pandas prefered format
        p_df = p_df.set_index('datetime').sort_index() # Index on time
        df_sub = p_df.loc[start:end] # Get the subset



    #Spatial Subset
    if bounding != None:
        #if both a scan_time and bounding area is given then the spatial subset start from the already
            #temporally subset dataframe... if not will start with the full platform dataframe
        if scan_time != None:
            aaa= df_sub.loc[(p_df['lat']>=bounding['ymin'])]
        else:
            aaa = p_df.loc[(p_df['lat']>=bounding['ymin'])]
        bbb = aaa.loc[(aaa['lat']<=bounding['ymax'])]
        ccc = bbb.loc[(bbb['lon']>=bounding['xmin'])]
        df_sub = ccc.loc[(ccc['lon']<=bounding['xmax'])]
        if print_long == True: print('Dataset has been spatially subset')

    #Test to ensure that there is valid data in the subrange
        #(aka the platform was deployed during the time of radarscan or that any of the points fall within the map area)
#    try:
#        p_test= df_sub.iloc[1]
#        p_deploy = True
#    except:
#        p_deploy = False
#        error_printing(e_test)

    # test to ensure that there is valid data in the subrange (aka the platform was deployed during the time of radarscan)
    p_deploy = (df_sub.size > 0)

    if not p_deploy:
        print('The platform was not deployed at this time')

    return df_sub, p_deploy

# Cache backed version
cached_grab_platform_subset = function_cache_memory.cache( grab_platform_subset )


# * * * * * *
def platform_plot(file, radartime, var, ax_n, p_attr, print_long, e_test, border_c='xkcd:light grey',labelbias=(0,0)):
    ''' Plot the in situ platform markers, barbs and pathline
    ----
    INPUTS:
    file: the in situ pandas dataframe
    var: dictionary containing info relatinto p_var (ie name, max, min)
    p_attr: dictionary containing platform styling info (color, shape etc)
    radartime, ax: time of the radar scan, axes of the subplot to plot on

    Optional Inputs:
    border_c, labelbias: color of background for the pathline, if you want to add labels directly to the plot this can offset it from the point
    '''


    #grab the subset of data of +- interval around radar scan
    p_sub, p_deploy = cached_grab_platform_subset(file, print_long, e_test, scan_time=radartime)

    if p_deploy == True:
        #plot the line that extends +/- min from the platform location
        print("p_sub = ", p_sub)
        plot_colourline(p_sub['lon'],p_sub['lat'],p_sub[var['name']],cmocean.cm.curl,ax_n,ccrs.PlateCarree(),color=border_c,amin=var['min'],amax=['max'])

        #plot the platform marker at the time closest to the scantime (aka the time at the halfway point of the subset platform dataframe)
        ax_n.plot(p_sub.iloc[[int(len(p_sub.index)/2)],['lon']], p_sub.iloc[[int(len(p_sub.index)/2)], ['lat']],ttransform=ccrs.PlateCarree(), marker=p_attr['m_style'],markersize=20, markeredgewidth='3',color=p_attr['m_color'],path_effects=[patheffects.withStroke(linewidth=12,foreground='k')], zorder=10)
        #plot labels for the marker on the plot itself
        #d2=plt.text(lon_sub[int(len(lon_sub)/2)]+labelbias[0], lat_sub[int(len(lat_sub)/2)]+labelbias[1], file[58:62],#transform=ccrs.PlateCarree(), fontsize=20, zorder=9, path_effects=[patheffects.withstroke(linewidth=4,foreground=color)])

        #plot a dot at the end of the colorline in the direction the platform is moving (aka the last time in the subset dataframe)
        ax_n.plot(p_sub.iloc[[-1],['lon']], p_sub.iloc[[-1],['lat']], transform=ccrs.PlateCarree(), marker='.',markersize=10, markeredgewidth='3',color='k',zorder=9)

        #plot windbarbs
        try:
            #p_sub[::x,'col_name'] returns every x'th value of th column data
            stationplot = metpy.plots.StationPlot(ax_n, p_sub.iloc[::30,'lon'], p_sub.iloc[::30,'lat'], clip_on=True, transform=ccrs.PlateCarree())
            stationplot.plot_barb(p_sub.iloc[::30,'U'], p_sub.iloc[::30,'V'], length=7)
        except: error_printing(e_test)

    elif p_deploy == False:
        if print_long==True: print('The platform was not deployed at this time')

    if print_long==True: print('made it through platform_plot')
    return



##  #  #  #  #  #  #  #  #  #  #  #  #  # # # #  #  #  #  #  #  #  #  #  #  # #  #  #  # # #  #  #  #  #  #  #  #  # #  #
##  #  #  #  #  #  #  #  #  #  #  #  #  # # # #  #  #  #  #  #  #  #  #  #  # #  #  #  # # #  #  #  #  #  #  #  #  # #  #
def ppiplot(r_only,radar,download_dir,day,globalamin, globalamax,p_var,p_of_int,CS3,e_test):
    if print_long== True: print('~~~~~~~~~~~Made it into ppiplot~~~~~~~~~~~~~~~~~~~~~')

    SMALL_SIZE, MS_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 23, 28, 33, 50
    plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=MS_SIZE)       # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    ## Determine Sweep for radar data
    swp=0

    ## Gather radar metadata
    # Find the beginning location of designated sweep in data array
    index_at_start = radar.sweep_start_ray_index['data'][swp]
    # Time at which scanning began
    currentscantime= num2date(radar.time['data'][index_at_start],radar.time['units'])
    print(currentscantime)
    # Convert time into fancy date string to use in title
    fancy_date_string_utc = currentscantime.strftime('%Y-%m-%d %H:%M UTC')
    print(fancy_date_string_utc)

    ## Determine elevation tilt for a chosen sweep and round to a decimal place
    elev=np.round(radar.get_elevation(swp)[0],1)

    ## Determine the 4letter radar ID
    r_site=radar.metadata['instrument_name']
    print('Gathered Radar Data')

    ## Set up plot size
    if r_only == True:
        fig = plt.figure(figsize=(10,8),facecolor='white')
        gs= gridspec.GridSpec(nrows=2,ncols=1)
    else:
        ## Set up figure
        fig = plt.figure(figsize=(32,20),facecolor='white')

        ## Establish the gridspec layout
        gs= gridspec.GridSpec(nrows=2, ncols=5,width_ratios=[.5,8,1,8,.25],height_ratios=[3,2],wspace=.25,hspace=.1)

        ## There are extra columns which are spacers to allow the formating to work
        # Uncomment to visualize the spacers
        #  ax_c=fig.add_subplot(gs[0,2])
        #  ax_l=fig.add_subplot(gs[0,0])
        #  ax_s=fig.add_subplot(gs[0,4])
    display = pyart.graph.RadarMapDisplay(radar)

    ## Plot instuments on subplots
    post=radar_subplots(radar,'refl',swp,fig,display,currentscantime,globalamin,globalamax,p_var,p_of_int,gs,r_site,e_test)
    post=radar_subplots(radar,'vel',swp,fig,display,currentscantime,globalamin,globalamax,p_var,p_of_int,gs,r_site,e_test)
    if r_only == False:
        time_series(day,fig,currentscantime,globalamin, globalamax, p_var,gs,radar_sub=True)

    ## Plot platform colorbar
    if p_var == "Thetae":
        c_lab= "Equivalent Potential Temp [K]"
    elif p_var == "Thetav":
        c_lab= "Virtual Potential Temp [K]"
    #print(post)
    cbar_ax = plt.axes([.514, post.y0,.014, post.y1-post.y0])#left, bottom, width, height
    cbar=plt.colorbar(CS3,cax=cbar_ax, orientation='vertical', label=c_lab, ticks=MaxNLocator(integer=True))#,ticks=np.arange(globalamin,globalamax+1,2))

    ## Plot title
    title=r_site+' '+str(elev)+r'$^{\circ}$ PPI '+ fancy_date_string_utc
    plt.suptitle(title,y=.92)

    save_dir = config.g_plots_directory + day + '/radar/Nexrad/plots/'

    if not os.path.exists(save_dir):
        Path(save_dir).mkdir(parents=True, exist_ok=True)


    ## Finish plot
    plt.savefig(save_dir+currentscantime.strftime('%m%d%H%M')+'_'+r_site+'_'+p_var+'.png' ,bbox_inches='tight',pad_inches=.3)
    plt.close()

    if print_long== True: print('~~~~~~~~~~~made it through ppiplot~~~~~~~~~~~~~~~~~~')
    print('Done Plotting')
    print('******************************************************************************************************')

    ## Makes a ding noise
    print('\a')
    return

def radar_subplots(radar, mom,swp,fig,display,currentscantime,globalamin,globalamax,p_var,p_of_int,gs,r_site,e_test):
    if print_long== True: print('~~~~~~~~~~~made it into radar_subplots~~~~~~~~~~~~~~')

    ## SET UP PLOTTING CONTROLS
    NSSLmm, NEBmm, UASd = True, True, False #which other platforms do you wish to plot
    KA_m, NOXP_m, ASOS_AWOS_m, MESONETS_m, METAR_m, WSR88D_m= True, False, False, True, False, True #additional markers to locate on the map
    country_roads, hwys, county_lines, state_lines = False, False, False, False #background plot features
    legend_elements=[]
    ###########################

    ## SET UP VARS FOR EACH RADAR MOMENTS
    if mom == 'refl':
        row, col, swp_id, feild = 0, 1, swp, 'reflectivity'
        c_scale, c_label = 'pyart_HomeyerRainbow', 'Radar Reflectivity [dbz]'
        vminb, vmaxb = -10., 75.
        p_title, leg = 'Reflectivity', True
    elif mom == 'vel':
        row, col, swp_id, feild = 0, 3, swp+1, 'velocity' #velocity data contained in the sweep greater than the refl at same level
        c_scale, c_label = 'pyart_balance', 'Velocity [m/s]'
        vminb, vmaxb = -40., 40.
        p_title, leg = 'Radial Velocity', False
    else:
        print("hey what just happened!\n")
        exit

    ## Bounding box, x km away from the probe of interest in all directions
    ymin,ymax,xmin,xmax,bndry = getLocation(p_of_int,currentscantime, e_test, offsetkm=40)

    # SET UP SUBPLOTS
    ax_n = fig.add_subplot(gs[row,col], projection=display.grid_projection)
    ax_n.text(.5,-.065,p_title, transform = ax_n.transAxes, horizontalalignment = 'center', fontsize = 40) #the radar subplot titles
    #  display.plot_ppi_map(feild, swp_id, title_flag=False,colorbar_flag=False,cmap=c_scale,ax=ax_n,vmin=vminb,vmax=vmaxb,embelish=False)
    display.plot_ppi_map(feild, swp_id, title_flag=False, colorbar_flag=False, cmap=c_scale,
                         ax=ax_n,vmin=vminb,vmax=vmaxb,
                         min_lon=xmin,max_lon=xmax,min_lat=ymin,max_lat=ymax,
                         embelish=False)

    ## PLOT RADAR MARKERS
    if WSR88D_m == True:
        #add markers for the WSR88D site
        #m_style,m_color,l_color,leg_str,legend_elements=platform_attr('WSR88D',legend_elements,radar_m=True, r_s=True,rad_site=r_site)
        p_attr,legend_elements=platform_attr('WSR88D',legend_elements,radar_m=True, r_s=True)
        m_style = p_attr['m_style']
        m_color = p_attr['m_color']
        l_color = p_attr['l_color']
        leg_str = p_attr['leg_str']


        ax_n.plot(radar.longitude['data'],radar.latitude['data'],marker=m_style, transform=ccrs.PlateCarree(),color=m_color,markersize=18,markeredgewidth=5,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')],zorder=10)
    if KA_m == True:
        print('Ka:made it into the loop')
        #get radar location data
        try:
            ka1dep, klat1, klon1, head1, rhib1, rhie1= det_radar_deps(currentscantime,'ka1',e_test)
            print('made it through ka1 loc')
        except:
            print('did not make through ka1 loc')
            error_printing(e_test)
        try:
            ka2dep, klat2, klon2, head2, rhib2, rhie2= det_radar_deps(currentscantime,'ka2',e_test)
            print('made it through ka2 loc')
        except:
            print('did not make through ka2 loc')
            error_printing(e_test)

        #plot the markers
        try:
            if np.logical_and(klat2>ymin,np.logical_and(klat2<ymax,np.logical_and(klon2>xmin,klon2<xmax))):
                print('made it in into ka2 marker')
                p_attr,legend_elements=platform_attr('Ka2',legend_elements,radar_m=True, r_s=True)
                m_style = p_attr['m_style']
                m_color = p_attr['m_color']
                l_color = p_attr['l_color']
                leg_str = p_attr['leg_str']

                ax_n.plot(klon2,klat2,marker=m_style, transform=ccrs.PlateCarree(),color=m_color,markersize=18,markeredgewidth=5,path_effects=[PathEffects.withStroke(linewidth=15,foreground='k')],
                        zorder=10)
                #d2=ax_n.text(klon2+0.006, klat2-0.011,'Ka2',transform=ccrs.PlateCarree(), zorder=10, path_effects=[PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')]) #textlabel on plot
        except:
            print('did not make it into ka2 marker')
            p_attr,legend_elements=platform_attr('Ka2',legend_elements,radar_m=True, r_s=True)
            m_style = p_attr['m_style']
            m_color = p_attr['m_color']
            l_color = p_attr['l_color']
            leg_str = p_attr['leg_str']


            error_printing(e_test)
        try:
            if np.logical_and(klat1>ymin,np.logical_and(klat1<ymax,np.logical_and(klon1>xmin,klon1<xmax))):
                print('made it in into ka1 marker')
                p_attr,legend_elements=platform_attr('Ka1',legend_elements,radar_m=True, r_s=True)
                m_style = p_attr['m_style']
                m_color = p_attr['m_color']
                l_color = p_attr['l_color']
                leg_str = p_attr['leg_str']

                ax_n.plot(klon1,klat1,marker=m_style, transform=ccrs.PlateCarree(),color=m_color,markersize=18,markeredgewidth=5,path_effects=[PathEffects.withStroke(linewidth=15,foreground='k')],
                        zorder=10)
                #d1=ax_n.text(klon1+0.006, klat1+0.005, 'Ka1',transform=ccrs.PlateCarree(), zorder=10, path_effects=[PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')])
        except:
            print('did not make it into ka1 marker')
            p_attr, legend_elements=platform_attr('Ka1',legend_elements,radar_m=True, r_s=True)
            m_style = p_attr['m_style']
            m_color = p_attr['m_color']
            l_color = p_attr['l_color']
            leg_str = p_attr['leg_str']


            error_printing(e_test)
    if NOXP_m == True:
        print('To be filled in')

    ## PLOT OTHER MARKERS
    if MESONETS_m == True:
        print("To Be filled in ")
        #WTM=pd.read_csv(filesys+'radarproc/West_TX_mesonets.csv')
        mesonet_file = config.g_mesonet_directory + 'radarproc/West_TX_mesonets.csv'
        if os.path.exists(mesonet_file):
            WTM=pd.read_csv(mesonet_file)
            WTM.name="WTx Mesonet"
            #ax_n.plot(WTM['Lon'],WTM['Lat'], transform=ccrs.PlateCarree(), marker=r'$\AA$',markersize=10, color='k',zorder=9)
            p_attr,legend_elements=platform_attr(WTM.name,legend_elements,r_s=True)
            m_style = p_attr['m_style']
            m_color = p_attr['m_color']
            l_color = p_attr['l_color']
            leg_str = p_attr['leg_str']


            ax_n.plot('Lon','Lat',data=WTM, transform=ccrs.PlateCarree(), marker=m_style,linestyle='None', markersize=15, color='k')
            proj=display.grid_projection
            WTM_sub=grab_platform_area_subset(WTM,ymin,ymax,xmin,xmax,p_var)
            for x,y,lab in zip(WTM_sub['Lon'],WTM_sub['Lat'],WTM_sub['Stn_ID']):
                #  print(lab)
                #  print(x)
                #  print(y)
                #  print(lab)
                trans =proj.transform_point(x,y,ccrs.Geodetic())
                print(trans)
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
    if METAR_m == True:
        print("To Be filled in ")
    if ASOS_AWOS_m == True:
        print("To Be filled in ")


    radartime = currentscantime

    ## PLOT INSITU TORUS PLATFORMS
    if NSSLmm == True:
        for NSSLMM in NSSLMM_df:
            name = NSSLMM.get('name', 'fixmenow1')
            if print_long== True: print(name)
            P_Attr,legend_elements=platform_attr(name,legend_elements,print_long,r_s=True)

            var = {'name': 'fixthisname', 'max': 1000.0, 'min': 0.0}
            print("fix me I am broken A")
            #platform_plot(NSSLMM, radartime, var, ax_n, P_Attr, print_long, e_test )
    if NEBmm == True:
        for UNLMM in UNLMM_df:
            name = UNLMM.get('name', 'fixmenow1')
            if print_long== True: print(name)
            P_Attr,legend_elements=platform_attr(name,legend_elements,print_long,r_s=True)
            var = {'name': 'fixthisname', 'max': 1000.0, 'min': 0.0}
            print("fix me I am broken B")
            #platform_plot(UNLMM,radartime, var, ax_n, P_Attr, print_long,e_test)
    if UASd == True:
        for UAS in UAS_files:
            name = UAS.get('name', 'fixmenow1')
            P_Attr,legend_elements=platform_attr(name,legend_elements,print_long,r_s=True)
            var = {'name': 'fixthisname', 'max': 1000.0, 'min': 0.0}
            print("fix me I am broken C")
            #platform_plot(UAS,radartime, var, ax_n, P_Attr, print_long, e_test)
                         # border_c='xkcd:very pale green',labelbias=(0,0.01))


    ## DEAL WITH COLORBARS
    # Attach colorbar to each subplot
    divider = make_axes_locatable(plt.gca())
    c_ax = divider.append_axes("right","5%", pad="2%",axes_class=plt.Axes)
    sm = plt.cm.ScalarMappable(cmap=c_scale,norm=matplotlib.colors.Normalize(vmin=vminb,vmax=vmaxb))
    sm._A = []
    cb = plt.colorbar(sm,cax=c_ax,label=c_label)

    ## Get position information to pass along to the remaining colorbar
    post=ax_n.get_position()


    ## SET UP LEGENDS
    l_list=legend_elements.tolist()
    if leg == True: #this means you are currently making the left subplot
        #add legend for platform markers
        l=ax_n.legend(handles=l_list,loc='center right', bbox_transform=ax_n.transAxes, bbox_to_anchor=(0,.5), handlelength=.1,title="Platforms",shadow=True,fancybox=True,ncol=1,edgecolor='black')
        l.get_title().set_fontweight('bold')

    else: #this means you are currently making the right subplot
        pass

    if print_long== True: print('~~~~~~~~~~~Made it through radar_subplots~~~~~~~~~~~')
    return post

# * * * * * *
def time_series(day,fig,currentscantime,globalamin,globalamax,p_var,gs,radar_sub=False):
    if print_long== True: print('~~~~~~~~~~~Made it into time_series~~~~~~~~~~~~~~~~~')

    # Convert scan_time from netcdf cftime to python standard datetime
    currentscantime = pd.to_datetime(currentscantime.strftime())

    if radar_sub == True:
        ax_n=fig.add_subplot(gs[1,:])
    else: #will only plot the time series (not the radar subplots)
        fig, ax_n= plt.subplots()

    ## MAKE THE TIME SERIES
    for platform_df in [NSSLMM_df,UNLMM_df,UAS_df]:
        #can remove this if/else statement once you fill in the UAS portion
        if platform_df == UAS_df:
            # Need to come back and fill in UAS
            # Can remove this if/else statement once you fill in the UAS portion
            pass
        else:
            for platform_file in platform_df:
                name = platform_file.get('name', 'fixmenow2')
                if print_long== True: print(str(name))

                ## Set up line colors and legend labels
                b_array=[] #blank array
                p_attr, legend_elements = platform_attr(name,b_array) #get the attributes (color, label, etc)

                ## Should ploting data be masked?
                # if you only want to mask certain platforms set up an if-statement to assign the "masking" variable
                masking = True #True/False variable
                plot_data= maskdata(p_var,platform_file,masking)

                ## Plot
                ax_n.plot(platform_file['datetime'],plot_data,linewidth=3,color=p_attr['l_color'],label=p_attr['leg_str']) #assigning label= is what allows the legend to work

    ## Set up XY axes tick locations
    ax_n.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M')) #strip the date from the datetime object
    ax_n.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax_n.xaxis.set_minor_locator(AutoMinorLocator(6)) # set up minor ticks (should be a multiple of ten intervals (ie 10,20,30... min spans)
    ax_n.yaxis.set_major_locator(MultipleLocator(5)) # set up major tick marks (this is set up to go by 5's will want to change for diff vars)
    ax_n.yaxis.set_minor_locator(AutoMinorLocator(5)) # set up minor ticks (this have it increment by 1's will want to change for diff vars)

    ## Set up axes formats (lables/ ticks etc)
    if p_var == "Thetae":
        ylab = "Equivalent Potential Temp [K]"
    elif p_var == "Thetav":
        ylab = "Virtual Potential Temp [K]"
    ax_n.set_ylabel(ylab)
    ax_n.set_xlabel('Time (UTC)')
    ax_n.tick_params(which='major', width=2,length=20, color='black')
    ax_n.tick_params(which='minor', width=2,length=10, color='grey')

    ## Set up grid for plot
    ax_n.grid(True)
    ax_n.grid(which='major', axis='x', color='grey', linewidth=2.5)
    ax_n.grid(which='minor', axis='x')
    ax_n.grid(which='major', axis='y', linestyle='--', linewidth=2)
    ax_n.grid(which='minor', axis='y', linestyle=':')

    ## Set up vertical lines that indicate the time of radar scan and the timerange ploted (via filled colorline) on the radar plots
    if radar_sub == True:
        ax_n.axvline(currentscantime,color='r',linewidth=4, alpha=.5, zorder=1)
        ax_n.axvspan(currentscantime-timedelta(minutes=5),currentscantime+timedelta(minutes=5), facecolor='0.5', alpha=0.4)

    ## Include the legend
    leg=ax_n.legend()
    for line in leg.get_lines():
        line.set_linewidth(12)

    if print_long== True: print('~~~~~~~~~~~Made it through time_series~~~~~~~~~~~~~~')
    return

#*****************************
#****************************************

def daterange(start_time, duration, iteration):
    ''' Creates an array of dates given a start time and length of time
        start_time: datetime string
        duration:   timedelta for a given amount of time
        iteration:  increment between datetimes'''
    i = int(duration.total_seconds()/iteration.total_seconds())
    j = iteration.total_seconds()/60.
    datetimes = []
    for n in range(i+1):
        datetimes.append(start_time+timedelta(minutes=float(j*n)))
    return datetimes

def circle(geod, lon, lat, radius, n_samples=360):
    """
    Return the coordinates of a geodetic circle of a given
    radius about a lon/lat point.

    Radius is in meters in the geodetic's coordinate system.
    """
    lons, lats, back_azim = geod.fwd(np.repeat(lon, n_samples),      #lons
                                     np.repeat(lat, n_samples),      #lats
                                     np.linspace(360, 0, n_samples), #azimuth
                                     np.repeat(radius, n_samples),   #distance
                                     radians=False,
                                     )
    return lons, lats, back_azim

def rring(lon, lat, radius_km, display, rhis=False,
          rhi_angles=[], heading=0):
    geod = Geod(ellps='WGS84')
    n_samples = 100
    geoms = []
    # Find lats/lons for range circle
    lons, lats, azim = circle(geod, lon, lat, radius_km * 1e3, n_samples)
    # Use shapely to create polygon a.k.a circle defined by lats/lons
    geoms.append(sgeom.Polygon(zip(lons, lats)))
    plt.gca().add_geometries(geoms, ccrs.PlateCarree(), facecolor='none',
                             edgecolor='k', alpha=0.75, linewidth=1)
    # Plot RHI "spokes" around circle
    if rhis==True:
        # find RHI azimuths that are north-relative
        rhi_angles = (np.asarray(rhi_angles)+heading) % 360
        geoms2 = []
        # Plot a line "spoke" for each RHI azimuth listed
        for angle in rhi_angles:
            ang_loc = np.argmin(np.abs((azim+180)-angle))
            rhi_lon = lons[ang_loc]
            rhi_lat = lats[ang_loc]
            geoms2.append(sgeom.LineString([(rhi_lon,rhi_lat),(lon,lat)]))
            plt.gca().add_geometries(geoms2, ccrs.PlateCarree(),
                                     edgecolor='k', alpha=0.75, linewidth=1)
    return

def plot_nexrad_loc(station,start_time,end_time,kloc_file,dir_out,radar_1=True,radar_2=False):
    ''' Plots NEXRAD reflectivity and velocity images relative to the
        location of the TTU Ka-band mobile radar.

        INPUT
        station:    4-letter WSR-88D station name (string)
        start_time: time at which to begin radar images (datetime object)
        end_time:   time at which to end radar images (datetime object)
        kloc:       filename of csv used to read in radar locations/times
                    format:'time_begin'; beginning of deployment
                           'time_end'; end of deployment
                           'lat','lon'; loction of deployment
                           'heading'; direction (in degrees) the radar was pointed
                           'radar'; "ka1" or "ka2"
                           'rhi_angles'; list of azimuths from deployment
    '''
    #print(kloc_file)
    if radar_1==True:
        ka1_loc_file = kloc_file[0]
        ka1_loc = pd.read_csv(ka1_loc_file,
                            parse_dates=['time_begin'])
        ka1_loc['time_end']=pd.to_datetime(ka1_loc['time_end'])
    if radar_2==True:
        ka2_loc_file = kloc_file[1]
        ka2_loc = pd.read_csv(ka2_loc_file,
                            parse_dates=['time_begin'])
        ka2_loc['time_end']=pd.to_datetime(ka2_loc['time_end'])

    # Grab the radar data between start time and end time
    print(start_time, end_time)
    radar_namelist, radar_list = get_radar_from_aws(station, start_time, end_time)

    # Loop through each radar volume to plot
    for radar in radar_list:
        sweep=0
        # Find the beginning location of designated sweep in data array
        index_at_start = radar.sweep_start_ray_index['data'][sweep]
        # Time at which scanning began
        time_at_start_of_radar = num2date(radar.time['data'][index_at_start],
                                          radar.time['units'])
        # Convert time into fancy date string to use in title
        eastern = pytz.timezone('US/Central')
        local_time = eastern.fromutc(time_at_start_of_radar)
        fancy_date_string = local_time.strftime('%A %B %d at %I:%M %p %Z')
        fancy_date_string_utc = time_at_start_of_radar.strftime('%Y-%m-%d %H:%M UTC')

        print('Gathered Radar Data')

        # Grab the location data for both radars
        # ka1_loc = k_loc[k_loc['radar']=='ka1'].reset_index()
        # ka2_loc = k_loc[k_loc['radar']=='ka2'].reset_index()

        # Begin figure
        fig = plt.figure(figsize = [20,8])
        display = pyart.graph.RadarMapDisplay(radar)
        proj = display.grid_projection
        lat_0 = display.loc[0] # Use radar data extent to define plot region
        lon_0 = display.loc[1]
        lat_1 = 33.96

        # Find location and heading of Ka-1 during NEXRAD scan
        if radar_1 == True:
            depl_time1  = np.abs(ka1_loc['time_begin']-time_at_start_of_radar).argmin()
            start_time1 = ka1_loc['time_begin'].iloc[[depl_time1]]
            end_time1   = ka1_loc['time_end'].iloc[[depl_time1]]
            timestr     = np.str(time_at_start_of_radar)
            klat1,klon1 = ka1_loc[['lat','lon']].iloc[[depl_time1]].values[0]
            ka1_heading = ka1_loc['heading'][depl_time1]
            # print('test_1',depl_time1, start_time1, end_time1)
        if radar_2 == True:
            depl_time2  = np.abs(ka2_loc['time_begin']-time_at_start_of_radar).argmin()
            start_time2 = ka2_loc['time_begin'].iloc[[depl_time2]]
            end_time2   = ka2_loc['time_end'].iloc[[depl_time2]]
            klat2,klon2 = ka2_loc[['lat','lon']].iloc[[depl_time2]].values[0]
            ka2_heading = ka2_loc['heading'][depl_time2]
            # print('test_2',depl_time2, start_time2, end_time2)

        # NEED TO CHANGE THIS...RHI ANGLES ARE HARD CODED IN RIGHT NOW
        #rhib = ka2_loc['rhib'][depl_time1]
        #rhie = ka2_loc['rhie'][depl_time1]
        #rhi_angles = np.arange(rhib, rhie+10,10)
        # print('rhib',rhib,'rhie',rhie,'range',rhi_angles)

        # Plot reflectivity on the left
        ax1 = fig.add_subplot(1,2,1,projection=proj)
        plt.tight_layout()
        field='reflectivity'
        title = field+' \n'+fancy_date_string+' ('+fancy_date_string_utc+') '
        offset = 0.75
        ax1.plot(klon1,klat1,color='r',markersize=200,transform=proj)
        display.plot_ppi_map(
                field, 0, colorbar_flag=True,
                title=title, ax=ax1,
                vmin=-10, vmax=75,
                min_lon=klon1-offset,max_lon=klon1+offset,
                min_lat=klat1-offset,max_lat=klat1+offset)

        if radar_1==True:
            cond11 = (time_at_start_of_radar-start_time1)[depl_time1]>dt.timedelta(minutes=0)
            cond21 = ((end_time1-time_at_start_of_radar)[depl_time1] >dt.timedelta(minutes=0))
        if radar_2==True:
            cond12 = (time_at_start_of_radar-start_time2)[depl_time2]>dt.timedelta(minutes=0)
            cond22 = ((end_time2-time_at_start_of_radar)[depl_time2] >dt.timedelta(minutes=0))

        if radar_1==True:
            if cond21:
                if cond11:
                    # Mark the location of Ka-1
                    display.plot_point(klon1,klat1,marker='+',color='k',
                                       label_text='Ka2',label_offset=(0.01,0.02),
                                       markersize=8,markeredgewidth=2)
                    # Plot a circle indicating the range of Ka1
                    rring(klon1,klat1, 15, display, rhis=False)
                          # rhi_angles=rhi_angles, heading=ka1_heading)
        if radar_2 == True:
            if cond22:
                if cond12:
                    # Mark the location of Ka-2
                    display.plot_point(klon2,klat2,marker='+',color='k',
                                       label_text='Ka1',label_offset=(0.01,-0.03),
                                       markersize=8,markeredgewidth=2)
                    # Plot a circle indicating the range of Ka2
                    rring(klon2,klat2, 15, display, rhis=False)
                          #rhi_angles=rhi_angles, heading=ka2_heading)

        # Plot velocity on the right
        ax2 = fig.add_subplot(1,2,2,projection=proj)
        plt.tight_layout()
        field='velocity'
        title = field+' \n'+fancy_date_string+' ('+fancy_date_string_utc+') '
        norm, cmap = ctables.registry.get_with_steps('NWSVelocity', 16, 16)
        display.plot_ppi_map(
                field, 1, colorbar_flag=True,
                title=title,#+obtitle,
                cmap=cmap, ax=ax2,
                vmin=-20, vmax=20,
                min_lon=klon1-offset,max_lon=klon1+offset,
                min_lat=klat1-offset,max_lat=klat1+offset)

        if (radar_1 == True) and (cond21):
            if cond11:
                # Mark the location of Ka-1
                display.plot_point(klon1,klat1,marker='+',color='k',
                                   label_text='Ka1',label_offset=(0.01,0.02),
                                   markersize=8,markeredgewidth=2)
                # Plot a circle indicating the range of Ka1
                rring(klon1,klat1, 15, display, rhis=False)
                      # rhi_angles=rhi_angles, heading=ka1_heading)
        if (radar_2 == True) and (cond22):
            if cond12:
                # Mark the location of Ka-2
                display.plot_point(klon2,klat2,marker='+',color='k',
                                   label_text='Ka2',label_offset=(0.01,-0.03),
                                   markersize=8,markeredgewidth=2)
                # Plot a circle indicating the range of Ka2
                rring(klon2,klat2, 15, display, rhis=False)
                      #rhi_angles=rhi_angles, heading=ka2_heading)

        #plt.savefig("test.png")
        title=(dir_out+'nexrad_'+timestr[0:4]+timestr[5:7]+timestr[8:10]+'_'+timestr[11:13]+timestr[14:16]+'.png')
        plt.savefig(str(title),format='png')
        #plt.savefig(dir_out+'nexrad_'+timestr[0:4]+timestr[5:7]+timestr[8:10]+'_'+timestr[11:13]+timestr[14:16]+'.png')
        print(timestr)
        plt.close()
    return

def radar_from_nextrad_file(radar_file):
    radar = pyart.io.read_nexrad_archive(radar_file)
    return radar

# Note: Cached version is cached on the file name, not the file contents.
# If file contents change you need to invalidate the cache or pass in the file contents directly to this function
cached_radar_from_nextrad_file = function_cache_memory.cache( radar_from_nextrad_file )


def plot_radar_file(radar_file, r_only, config, day, globalamin, globalamax, p_var, p_of_int,CS3, e_test):
    i = 0
    radar = None

    print("[{}] open_pyart, scan file_name = {}\n".format(i, radar_file))
    try:
        radar = cached_radar_from_nextrad_file (radar_file)
    except:
        print("****************    Failed to convert file: ", radar_file)
        print("****************    Failed to convert file: ", radar_file)
        print("****************    Failed to convert file: ", radar_file)

    #testing_plots(radar,i,j,r_only,globalamin,globalamax,p_var,e_test)
    #plot_with_WSR(radar,i,j)

    if radar is not None:
        ppiplot(r_only, radar, config.g_download_directory, day, globalamin,globalamax, p_var, p_of_int,CS3, e_test)


###############################################################################################
##########################################
# Create new datasets for each platform ##
##########################################
##Crop the mesonet timeframe?
#crop the start or end time to a time you specify (comment out to do the whole day)
tstart = dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),15,0,0)
tend = None

NSSLMM_df, UNLMM_df, max_array, min_array = [], [], [], []

for MM in ['FFld','LIDR','Prb1','Prb2','WinS','CoMeT1','CoMeT2','CoMeT3','UAS']:

    max_val = None
    min_val = None

    #load the NSSL mesonets
    if (MM in ['FFld','LIDR','Prb1','Prb2','WinS']):
        #try:
            MMfile= config.g_mesonet_directory + day +'/mesonets/NSSL/'+MM+'_'+day[2:]+'_QC_met.dat'

            files = glob.glob(MMfile)
            print("B: files: ", files)

            if len(files) > 0:
                file = files[0]  # Use first file only !!!!!! ??????

                df = read_nsslmm(file,tstart,tend)

                df.name = MM

                if MM == 'Prb1':
                    Prb1_df = df

                NSSLMM_df.append(df)

                masked_df = maskdata(p_var, df, mask=True)
                max_val, min_val= masked_df.max(), masked_df.min()

            if MM == probe_of_interest:
                p_of_int=NSSLMM_df[-1]

        #except: error_printing(e_test)

    #load the UNL MM files
    elif (MM in ['CoMeT1', 'CoMeT2', 'CoMeT3']):
        #try:
            MMfile=config.g_mesonet_directory + day +'/mesonets/UNL/UNL.'+MM+'.*'
            files = glob.glob(MMfile)
            print("C: files: ", files)

            if len(files) > 0:
                file = files[0]  # Use first file only !!

                df = read_unlmm(file,tstart,tend)

                df.name = MM

                UNLMM_df.append(df)

                masked_df = maskdata(p_var, df, mask=True)
                max_val, min_val= masked_df.max(), masked_df.min()

            if MM == probe_of_interest:
                p_of_int=UNLMM_df[-1]

        #except: error_printing(e_test)

    #load the UAS files
    elif (MM in ['UAS','UAS_fillers']):
        UAS_df='Will come back and fill in'


    if max_val is not None:
        max_array.append(max_val)

    if min_val is not None:
        min_array.append(min_val)

print("max_array= ", max_array)
print("min_array= ", min_array)

#Alternate global variable max/min for plotting purposes
globalamax = np.max(max_array)
globalamin = np.min(min_array)

#######################################################################################################
#Dummy plot for scatter
# Make a scale covering max and min
cmap=cmocean.cm.curl
Z = [[0,0],[0,0]]
levels = np.arange(globalamin,globalamax+1,1)
CS3 = plt.contourf(Z, levels, cmap=cmap)#,transform=datacrs)
plt.clf() # clear current figure
# ************************************


# ???? is Prb1_df really what you think it is?

#print(Prb1_df)
#for i in len(Prb1_df):
#    Prb1_df['Radar_ID']=nearest_WSR(
test=det_nearest_WSR(Prb1_df)
r_ofintrest=test.Radar_ID.unique()

timeranges_each_r=pd.DataFrame()
for r in r_ofintrest:
    print(r)
    t=test.loc[test.Radar_ID==r, ['datetime']].rename(columns={'datetime':r})
    ts,te=t.min(),t.max()
    timeranges_each_r=pd.concat([timeranges_each_r,t],axis=1)
    print("start ",ts[r])
    print("end ",te[r])
    print("***")

    radar_files = get_WSR_from_AWS(ts[r],te[r],r,config.g_download_directory)

    print("Radar files to process:")
    print(radar_files)


    Parallel(n_jobs=5, verbose=10)(delayed(plot_radar_file)(radar_file, r_only, config, day, globalamin, globalamax, p_var, p_of_int,CS3, e_test) for radar_file in radar_files)

    #open the downloaded files as pyart objects
    #radar_list=[]
    #for radar_file in radar_files:
    #    i = i+1
    # plot_radar_file(i, radar_file, r_only, config, day, globalmin, globalmax, p_var, p_of_int,CS3, e_test)



   #plot_with_WSR(test2)
    #print(timeranges_each_r)

    #tt=t.rename(columns={'datetime':r})
    #tt=t['datetime'].to_numpy()
    #print(tt)
    #print(np.shape(tt))
    #print(tt)
    #timeranges_each_r[r]=pd.Series(tt)
    #timeranges_each_r.loc[r]=tt
print(timeranges_each_r)
