#import needed modules
######################
import matplotlib
matplotlib.use('agg')
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import xarray as xr
from xarray.backends import NetCDF4DataStore
import matplotlib.pyplot as plt
import cmocean
import gc
import os.path
from datetime import datetime
from datetime import date
import datetime as dt
from matplotlib.colors import ListedColormap
import matplotlib.pylab as pl
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from datetime import timedelta
from mpl_toolkits import axes_grid1
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import ShapelyFeature,NaturalEarthFeature
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
from cartopy.io.shapereader import Reader
import matplotlib.patheffects as PathEffects
import glob
import pyart
import os
import cartopy.io.shapereader as shpreader
from matplotlib.colors import Normalize
from metpy.plots import ctables
import metpy
import metpy.calc as mpcalc
from metpy.units import units
from scipy import ndimage
import pandas as pd
import osmnx as ox
import shutil
from os import path
from scipy import interpolate
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.lines import Line2D
import read_platforms
import sys, traceback
from collections import namedtuple

##########################################################################
## VARIABLES
############
day = '20190524' #'YYYYMMDD'
radar_plotting= True #would you like to plot radar images? (set to False to only plot timeseries)
r_only=False #set to true if you want the radar plots only 
p_var = "Thetav" #which var to plot (current options; Thetae, Thetav)
filesys='/Users/severe2/Research/'
temploc='/Volumes/Samsung_T5/Research/TORUS_Data/'

## Read in defns that I have stored in another file (for ease of use/consistancy accross multiple scripts
from shared_defns import read_platforms, maskdata, platform_attr, error_printing

#Troubleshooting y/n
####################
#Set to True to print statements when entering/leaving definitions (helpful to understanding whats going on in the code but prints alot of lines to your terminal)
print_long = True# True/False variable

#Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting 
    ## To do do calls on defn from the shared_defns folder
    ## there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = True#True/False variable

##########################################################################
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
    
    for i in range(0,len(range_mask[:,0])):
        range_mask[i,:] = radar.range['data']>(radar.range['data'][-1]-1000.)
    
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

    return
    
# * * * * * * * 
def det_radar_deps(currentscantime, r_testing, e_test):
    ''' Determine if a given (ka) radar is deployed and if so assign the correct values to it. 
        r_testing is the name of the radar you are testing to see if deployed
    '''
    #test if radar deployed
    try:
        kadep=pd.read_csv(filesys+'TORUS_Data/'+day+'/radar/TTUKa/csv/'+day+'_deployments_'+r_testing+'.csv')
        
        #get radar positioning data
        k_lat, k_lon, head_ing, rhi_b, rhi_e = getRadarPositionData(currentscantime, kadep)
    
    except: 
        error_printing(e_test)
        #If one of the radars did not deploy at all on a given day this will set its values to np.nan
        k_lat, k_lon, head_ing, rhi_b, rhi_e = np.nan, np.nan, np.nan, np.nan, np.nan
    
    return k_lat, k_lon, head_ing, rhi_b, rhi_e

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


###########################
## PLOTTING DEFINTIONS  ###
###########################
def getLocation(current_lat, current_lon, offsetkm, given_bearing= False):
    ''' 
    This definition has two functions:
        1) If no bearing is specified it will return a namedtuple containing the max/min lat/lons 
            to form a square surrounding the point indicated by lat1,lon1 by x km. 
        2) If a bearing is given then the defintion will return one set of lat/lon values 
             (end_lat and end_lon) which is the location x km away from the point lat1, lon1
             if an observer moved in a straight line the direction the bearing indicated.  
    ----
    INPUTS
    current_lon & current_lat: the starting position 
    offsetkm: distance traveled from the starting point to the new locations 
    given_bearing: True/False, are you given a specified direction of "travel"
    '''
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
        
        #set up a namedtuple object to hold the new info
        box_extent = namedtuple('box_extent', ['ymin','ymax','xmin','xmax'])
        box=box_extent(ymin=min_lat, ymax=max_lat, xmin=min_lon, xmax= max_lon)
        
        return box
    
    else:
        #if a bearing is provided 
        bearing = (given_bearing/ 90.)* np.pi / 2.

        new_lat = np.arcsin(np.sin(lat1) * np.cos(offsetkm/R) + np.cos(lat1) * np.sin(offsetkm/R) * np.cos(bearing))
        new_lon = lon1 + np.arctan2(np.sin(bearing)*np.sin(offsetkm/R)*np.cos(lat1),np.cos(offsetkm/R)-np.sin(lat1)*np.sin(new_lat))
        end_lon = 180.0 * new_lon / np.pi
        end_lat = 180.0 * new_lat / np.pi
        
        return end_lat, end_lon 

#* * * * * 
def createCircleAroundWithRadius(clat, clon, radiuskm,sectorstart,sectorfinish,heading):
    #    ring = ogr.Geometry(ogr.wkbLinearRing)
    latArray = []
    lonArray = []
    for bearing in range(int(heading+sectorstart),int(heading+sectorfinish)): #degrees of sector
        lat2, lon2 = getLocation(clat,clon,radiuskm,given_bearing=bearing)
        latArray.append(lat2)
        lonArray.append(lon2)
    
    return lonArray,latArray

# * * * * * * 
def rhi_spokes_rings(ka_info,display,print_long,e_test):
    ''' Plot the RHI spoke and ring for the KA radars
    '''
    if print_long== True: print('made it into rhi_spokes_rings')

    for ka in ['KA1','KA2']:
        try:
            #produce spoke and ring
            for j in range(int(ka_info[ka]['rhib']),int(ka_info[ka]['rhie'])+1,10):
                #  print(j)
                ang= ka_info[ka]['head'] + j
                #  print(ang)
                if ang > 360.:
                    ang=int(ang-360.)

                #this plots a circle that connects the spokes
                A,B = createCircleAroundWithRadius(ka_info[ka]['lat'],ka_info[ka]['lon'],(ka_info['Rfile'].range['data'][-1]-500.)/1000., ka_info[ka]['rhib'],ka_info[ka]['rhie']+1,ka_info[ka]['head']) 
                C,D = getLocation(ka_info[ka]['lat'], ka_info[ka]['lon'], (ka_info['Rfile'].range['data'][-1]-500.)/1000., given_bearing=ang)
                display.plot_line_geo(A,B,marker=None,color='grey',linewidth=.25) #this plots a circle that connects the spokes
                display.plot_line_geo([ka_info[ka]['lon'],D], [ka_info[ka]['lat'],C], marker=None, color='k', linewidth=0.5, linestyle=":")
                #if np.logical_and(C>ymin,np.logical_and(C<ymax,np.logical_and(D>xmin,D<xmax))):
                    #d1=plt.text(D, C, str(ang),horizontalalignment='center',transform=ccrs.PlateCarree(),fontsize=10,zorder=9,path_effects=([PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')])
        except: error_printing(e_test)
    if print_long== True: print('made it through rhi_spokes_rings')
    return

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
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=color, linewidth=10.5, transform=datacrs, zorder=3)
    for i in np.arange(len(x)-1): #This is the colorramp colorline
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=c[i], linewidth=7.5, transform=datacrs, zorder=4)
    
    return

# * * * * * *
def grab_platform_subset(p_df, print_long, e_test, bounding=None, scan_time= None, time_offset= 300):
    ''' This def will take a given dataset and subset it either spatially or temporially
            1) If scan_time is given the data will be subset temporally 
                    Grabs the observed thermo data +/- x seconds around radar scan_time.
                    The length of the subset is det by time_offset which is defaults to +/- 300 secs. 
                    This means that unless a different time_offeset is specified the subset will span 10 min. 
                **NOTE: This does not automatically control the gray vertical box on the timeseries**
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
        aaa = p_df.loc[(p_df['datetime']>=scan_time-dt.timedelta(seconds=time_offset))]
        df_sub = aaa.loc[(aaa['datetime']<=scan_time+dt.timedelta(seconds=time_offset))]
        if print_long == True: print('Dataset has been temporally subset')

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
    try:
        p_test= df_sub.iloc[1]
        p_deploy = True
    except:
        p_deploy = False
        error_printing(e_test)

    return df_sub, p_deploy

# * * * * ** *    
def platform_plot(file, radartime, var, p_attr, ax, print_long, e_test, border_c='xkcd:light grey',labelbias=(0,0)):
    ''' Plot the in situ platform markers, barbs and pathline
    ----
    INPUTS:
    file: the in situ pandas dataframe
    var: dictionary containing info relating to p_var (ie name, max, min)
    p_attr: dictionary containing platform styling info (color, shape etc)
    radartime, ax: time of the radar scan, axes of the subplot to plot on

    Optional Inputs:
    border_c, labelbias: color of background for the pathline, if you want to add labels directly to the plot this can offset it from the point  
    '''
    if print_long== True: print('made it into platform_plot')

    #grab the subset of data of +- interval around radar scan
    p_sub, p_deploy = grab_platform_subset(file, print_long, e_test, scan_time=radartime)

    if p_deploy == True:
        #plot the line that extends +/- min from the platform location 
        plot_colourline(p_sub['lon'].values, p_sub['lat'].values, p_sub[var['name']].values, cmocean.cm.curl, ax, ccrs.PlateCarree(), color=border_c, amin=var['min'], amax=var['max'])
       
        #find the value of the index that is halfway through the dataset (this will be the index associated with radar_scantime) 
        mid_point=(p_sub.index[-1]-p_sub.index[0])/2
        col_lon, col_lat =p_sub.columns.get_loc('lon'), p_sub.columns.get_loc('lat')
        col_U, col_V =p_sub.columns.get_loc('U'), p_sub.columns.get_loc('V')

        #plot the platform marker at the time closest to the scantime (aka the time at the halfway point of the subset platform dataframe)
        ax.plot(p_sub.iloc[mid_point,col_lon], p_sub.iloc[mid_point,col_lat], transform=ccrs.PlateCarree(), marker=p_attr['m_style'], markersize=20, markeredgewidth='3',color=p_attr['m_color'], path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')], zorder=10)
        #plot labels for the marker on the plot itself
        #  d2=plt.text(p_sub.iloc[mid_point,col_lon]+labelbias[0], p_sub.iloc[mid_point,col_lat]+labelbias[1], file[58:62],transform=ccrs.PlateCarree(), fontsize=20, zorder=9, path_effects=[patheffects.withstroke(linewidth=4,foreground=color)])
  
        #plot a dot at the end of the colorline in the direction the platform is moving (aka the last time in the subset dataframe) 
        ax.plot(p_sub.iloc[-1,col_lon], p_sub.iloc[-1,col_lat],transform=ccrs.PlateCarree(), marker='.',markersize=10, markeredgewidth='3',color='k',zorder=9)
  
        #plot windbarbs
        try:
            #p_sub.iloc[::x,col_index] returns every x'th value 
            stationplot = metpy.plots.StationPlot(ax, p_sub.iloc[::30,col_lon], p_sub.iloc[::30,col_lat], clip_on=True, transform=ccrs.PlateCarree())
            stationplot.plot_barb(p_sub.iloc[::30,col_U], p_sub.iloc[::30,col_V],sizes=dict(emptybarb= 0),length=7)
            pass
        except: error_printing(e_test)

    elif p_deploy == False: 
        if print_long==True: print('The platform was not deployed at this time')

    if print_long== True: print('made it through platform_plot')
    return



##  #  #  #  #  #  #  #  #  #  #  #  #  # # # #  #  #  #  #  #  #  #  #  #  # #  #  #  # # #  #  #  #  #  #  #  #  # #  #
##  #  #  #  #  #  #  #  #  #  #  #  #  # # # #  #  #  #  #  #  #  #  #  #  # #  #  #  # # #  #  #  #  #  #  #  #  # #  #
def ppiplot(ka_info, var, pform, r_only, filesys, day, print_long, e_test):
    ''' 
    Initial plotting defenition: sets up fig size, layout, font size etc and will call timeseries and radar subplots
    ----------
    INPUTS
    ka_info : dictionary
        contains info about the ka radars (ie time, locations, sweep, azimuth, and RHI angles)
    var: dictionary 
        contains info about the variable being overlayed (ie Thetae etc). Contains name and the max and min value 
    pform: dictionary
        contains the pandas dataframes of the initu platforms
    r_only : string
        true/false variable that indicates whether only radar should be plotted (aka no Timeseries if True)
        (I have not actually tried this yet .... but in theroy this should work (you would need to play with formating))
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
    if r_only == True:
        fig = plt.figure(figsize=(20,10), facecolor='white')
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
    display = pyart.graph.RadarMapDisplay(ka_info['Rfile'])

    ## Plot instuments on subplots
    post=radar_subplots('refl', fig, gs, display, ka_info, var, pform, filesys, print_long, e_test)
    post=radar_subplots('vel', fig, gs, display, ka_info, var, pform, filesys, print_long, e_test)
    
    ## Plot timeseries 
    if r_only == False:
        time_series(fig, gs, ka_info, var, pform, print_long, e_test)

    ## Plot platform colorbar
    if var['name'] == "Thetae":
        c_lab= "Equivalent Potential Temp [K]"
    elif var['name'] == "Thetav":
        c_lab= "Virtual Potential Temp [K]"
    cbar_ax = plt.axes([.514, post.y0,.014, post.y1-post.y0])#left, bottom, width, height
    cbar=plt.colorbar(CS3,cax=cbar_ax, orientation='vertical', label=c_lab, ticks=MaxNLocator(integer=True))#,ticks=np.arange(globalamin,globalamax+1,2))

    ## Plot title
    title = thefile.split("/")[-1][10:13] + ' '+str(np.around(ka_info['Azimuth'],2))+r'$^{\circ}$ PPI ' +ka_info['Time'].strftime("%m/%d/%y %H:%M")+ ' UTC'
    plt.suptitle(title,y=.92)
    
    ## Finish plot 
    plt.savefig(filesys+'TORUS_Data/'+day+'/mesonets/plots/'+thefile.split("/")[-1][10:13]+'_'+ka_info['Time'].strftime('%m%d%H%M')+'_'+var['name']+'.png' ,bbox_inches='tight',pad_inches=.3)
    print(filesys+'TORUS_Data/'+day+'/mesonets/plots/'+thefile.split("/")[-1][10:13]+'_'+ka_info['Time'].strftime('%m%d%H%M')+'_'+var['name']+'.png')
    plt.close()

    if print_long== True: print('~~~~~~~~~~~made it through ppiplot~~~~~~~~~~~~~~~~~~')
    print('Done Plotting')
    print('******************************************************************************************************')

    ## Makes a ding noise 
    print('\a')
    return 

# * * * * * *  *
def radar_subplots(mom, fig, gs, display, ka_info, var, pform, filesys, print_long, e_test):
    ''' Plots each of the radar subplots and calls for the markers to be plotted for all the additional platforms
    ----
    INPUTS
    mom: str, indicates which radar moment you would like to plot (ie Reflectivity, Velocity, Spectrum Width etc)
    fig, gs, & display: information about the plot 
    ka_info, var, & pform: dictionarys, as described in the ppiplot comment 
            **** NOTE: I need to come back and do more with pform****
    '''
    if print_long== True: print('~~~~~~~~~~~made it into radar_subplots~~~~~~~~~~~~~~')
    
    ## SET UP PLOTTING CONTROLS
    NSSLmm, NEBmm, UASd = True, True, False #which other platforms do you wish to plot 
    rhi_ring= True #do you want the rhi spokes
    country_roads, hwys, county_lines, state_lines = False, False, False, False #background plot features
    legend_elements=[]
    ###########################

    ## SET UP VARS FOR EACH RADAR MOMENTS
    if mom == 'refl':
        row, col, field= 0, 1, 'refl_fix'
        c_scale, c_label='pyart_HomeyerRainbow','Radar Reflectivity [dbz]'
        vminb, vmaxb= -30.,30.
        p_title, leg ='Reflectivity', True
    elif mom == 'vel':
        row, col, field= 0, 3, 'vel_fix'
        c_scale, c_label='pyart_balance','Velocity [m/s]'
        vminb, vmaxb= -40.,40.
        p_title, leg ='Radial Velocity', False
    else:
        print("hey what just happened!\n")
        exit
   
    ## Read in the primary radar variable
    Main_R=ka_info['Main_Radar']
   
    ## Bounding box, x km away from the radar in all directions 
    box=getLocation(ka_info[Main_R]['lat'],ka_info[Main_R]['lon'],offsetkm=21)    
    
    ## SET UP SUBPLOTS   
    ax_n=fig.add_subplot(gs[row,col],projection=display.grid_projection)
    ax_n.plot(ka_info[Main_R]['lon'],ka_info[Main_R]['lat'],transform=display.grid_projection)
    ax_n.text(.5,-.065,p_title,transform=ax_n.transAxes,horizontalalignment='center',fontsize=40) #the radar subplot titles
    display.plot_ppi_map(field, ka_info['Swp_ID'],title_flag=False,colorbar_flag=False, cmap=c_scale, ax=ax_n, vmin=vminb, vmax=vmaxb, min_lon=box.xmin, max_lon=box.xmax, min_lat=box.ymin, max_lat=box.ymax, embelish=False) 
        
    ## PLOT RHI SPOKES
    if rhi_ring == True:
        rhi_spokes_rings(ka_info, display, print_long, e_test) #plot the spokes for radars

    ## PLOT RADAR MARKERS
    for ka in ['KA1', 'KA2']:
        try:
            if np.logical_and(ka_info[ka]['lat']>box.ymin,np.logical_and(ka_info[ka]['lat']<box.ymax,np.logical_and(ka_info[ka]['lon']>box.xmin,ka_info[ka]['lon']<box.xmax))):
                P_Attr,legend_elements=platform_attr(ka,print_long,r_s=True,radar_m=True,l_array=legend_elements)

                ax_n.plot(ka_info[ka]['lon'],ka_info[ka]['lat'],marker=P_Attr['m_style'], transform=ccrs.PlateCarree(),color=P_Attr['m_color'],markersize=18,markeredgewidth=5,path_effects=[PathEffects.withStroke(linewidth=15,foreground='k')],
                        zorder=10)
                #d2=ax_n.text(ka_info[ka]['lon']+0.006, ka_info[ka]['lat']-0.011,'Ka2',transform=ccrs.PlateCarree(), zorder=10, path_effects=[PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')]) #textlabel on plot 
        except: error_printing(e_test)
    
    ## PLOT INSITU TORUS PLATFORMS              
    if NSSLmm == True:
        for NSSLMM in NSSLMM_df:
            if print_long== True: print(NSSLMM.name)
            P_Attr,legend_elements=platform_attr(NSSLMM, print_long, r_s=True, l_array=legend_elements)
            platform_plot(NSSLMM, ka_info['Time'], var, P_Attr, ax_n, print_long, e_test)
    if NEBmm == True:
        for UNLMM in UNLMM_df:
            if print_long== True: print(UNLMM.name)
            P_Attr,legend_elements=platform_attr(UNLMM, print_long, r_s=True, l_array=legend_elements)
            platform_plot(UNLMM, ka_info['Time'], var, P_Attr, ax_n, print_long, e_test)
    if UASd == True:
        for UAS in UAS_files:
            P_Attr,legend_elements=platform_attr(UAS, print_long, r_s=True, l_array=legend_elements)
            platform_plot(UAS, ka_info['Time'], var, P_Attr, ax_n, border_c='xkcd:very pale green', labelbias=(0,0.01))
 
    
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

    ## PLOT BACKGROUND FEATURES
    if country_roads == True:
        ox.config(log_file=True, log_console=True, use_cache=True)
        G = ox.graph_from_bbox(box.ymax,box.ymin,box.xmax,box.xmin)
        ox.save_load.save_graph_shapefile(G, filename='tmp'+str(0), folder=filesys+'radarproc/roads/', encoding='utf-8')
        fname = filesys+'radarproc/roads/tmp'+str(0)+'/edges/edges.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree(), edgecolor='gray', linewidth=0.5)
        ax_n.add_feature(shape_feature, facecolor='none')
        shutil.rmtree(filesys+'radarproc/roads/tmp'+str(0)+'/')
    if hwys == True:
        fname = filesys+'radarproc/roads/GPhighways.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree(), edgecolor='grey')#edgecolor='black')
        ax_n.add_feature(shape_feature, facecolor='none')
    if county_lines == True:
        fname = filesys+'radarproc/roads/cb_2017_us_county_5m.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree(), edgecolor='gray')
        ax_n.add_feature(shape_feature, facecolor='none', linewidth=1.5, linestyle="--")
    if state_lines == True:
        states_provinces = cartopy.feature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
        ax_n.add_feature(states_provinces, edgecolor='black', linewidth=2)
     
    if print_long== True: print('~~~~~~~~~~~Made it through radar_subplots~~~~~~~~~~~')
    return post

# * * * * * * * 
def time_series(fig, gs, ka_info, var, pform, print_long, e_test, radar_sub=True):
    ''' Plot the time series of p_var from the various instrument platforms 
    ----
    INPUTS
    fig & gs: relating to the plot layout 
    ka_info, var, & pform: dictionarys, as described in the ppiplot comment 
            **** NOTE: I need to come back and do more with pform****
    radar_sub: True/False, If False will only plot timeseries (no radar).. In theroy have not actually done this yet
    '''
    if print_long== True: print('~~~~~~~~~~~Made it into time_series~~~~~~~~~~~~~~~~~')

    if radar_sub == True:
        ax_n=fig.add_subplot(gs[1,:])
    
    else: #will only plot the time series (not the radar subplots)
        fig, ax_n= plt.subplots()
    
    ## MAKE THE TIMESERIES 
    for platform_df in [NSSLMM_df,UNLMM_df,UAS_df]: 
        #can remove this if/else statement once you fill in the UAS portion
        if platform_df == UAS_df:
            # Need to come back and fill in UAS
            # Can remove this if/else statement once you fill in the UAS portion
            pass
        else:
            for platform_file in platform_df:
                if print_long== True: print(str(platform_file.name))
                
                ## Set up line colors and legend labels 
                P_Attr = platform_attr(platform_file, print_long)

                ## Should ploting data be masked?
                # if you only want to mask certain platforms set up an if-statement to assign the "masking" variable
                masking = True #True/False variable
                plot_data= maskdata(var['name'],platform_file,masking)

                ## Plot
                ax_n.plot(platform_file['datetime'],plot_data,linewidth=3,color=P_Attr['l_color'],label=P_Attr['leg_str']) #assigning label= is what allows the legend to work  
                    
    ## Set up XY axes tick locations
    ax_n.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M')) #strip the date from the datetime object
    ax_n.xaxis.set_major_locator(mdates.AutoDateLocator()) 
    ax_n.xaxis.set_minor_locator(AutoMinorLocator(6)) # set up minor ticks (should be a multiple of ten intervals (ie 10,20,30... min spans)
    ax_n.yaxis.set_major_locator(MultipleLocator(5)) # set up major tick marks (this is set up to go by 5's will want to change for diff vars)
    ax_n.yaxis.set_minor_locator(AutoMinorLocator(5)) # set up minor ticks (this have it increment by 1's will want to change for diff vars)
    
    ## Set up axes formats (lables/ ticks etc)
    if var['name'] == "Thetae":
        ylab = "Equivalent Potential Temp [K]"
    elif var['name'] == "Thetav":
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
        ax_n.axvline(ka_info['Time'],color='r',linewidth=4, alpha=.5, zorder=1)
        ax_n.axvspan(ka_info['Time']-timedelta(minutes=5),ka_info['Time']+timedelta(minutes=5), facecolor='0.5', alpha=0.4)
    
    ## Include the legend
    leg=ax_n.legend()
    for line in leg.get_lines():
        line.set_linewidth(12)
    
    if print_long== True: print('~~~~~~~~~~~Made it through time_series~~~~~~~~~~~~~~')
    return


###############################################################################################
##########################################
# Create new datasets for each platform ##
##########################################
##Crop the mesonet timeframe?
#crop the start or end time to a time you specify (comment out to do the whole day) 
tstart = dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),15,0,0)
tend = None

Platform_Types=['NSSL','UNL','UAS']
pform={
    'NSSL':{'FFld':{'data':'', 'name':'', 'masked_data':'', 'min':'', 'max':''},
            'LIDR':{},
            'Prb1':{},
            'Prb2':{},
            'WinS':{}
        },
    'UNL':{'CoMeT1':{}, 
           'CoMeT2':{}, 
           'CoMeT3':{}
        },
    'UAS':{'UASFillers':{}
        }
    }

class Pname:
    
    def __init__(self, Name, p_var, tstart, tend, day, filesys, print_long, e_test, mask=True):
        '''setting the ojects attr anything with self.something= can be called using AAA.something in your code after you have initilized the object
            this defn will only be run once per object
        '''
        #set up filepath for the platform
        if Name in ['FFld','LIDR','Prb1','Prb2','WinS']:
            File=filesys+'TORUS_data/'+day+'/mesonets/NSSL/'+Name+'_'+day[2:]+'_QC_met.dat'
            Ptype='NSSL'
        elif Name in ['CoMeT1', 'CoMeT2', 'CoMeT3']:
            File='/Volumes/Samsung_T5/Research/TORUS_Data/'+day+'/mesonets/UNL/UNL.'+Name+'.*'
            Ptype = 'UNL'
        elif Name in ['UAS']:
            print('Code for platform has not been written yet')
        else:
            print('There seems to have been error please check spelling and retry')
        
        #read in and intilize the dataset
        self.df= read_platforms(Ptype, File, tstart, tend)
        self.df.name='{}_df'.format(Name) #assign a name to the pandas dataframe itself 
        
        #get style info for the platforms marker
        self.m_style, self.m_color, self.l_color, self.leg_str, self.leg_entry= platform_attr(Name, print_long) 
        
        # apply mask (if applicable) and find the min and max of p_var for this platform 
        if mask== True: #mask is a bool, sent to True to by default
            self.mask_df= maskdata(p_var,self.df) #add masked dataset to the object 
            self.Min= self.mask_df.min()#use the masked dataset to find max and min of p_var
            self.Max= self.mask_df.max()
        elif mask == False:
            self.Min= self.df.min()[p_var] #use the orig dataset to find max/min
            self.Max= self.df.max()[p_var]


FFld= Pname('FFld',p_var,tstart,tend,day, filesys, print_long,e_test)    
print(FFld.Min)


#######
for p_type, p_info in pform.items():
    for pname in p_info:
        if print_long== True: print(pname)
        
        #Read in the data 
        if p_type == 'NSSL':
            try:
                File=filesys+'TORUS_data/'+day+'/mesonets/NSSL/'+pname+'_'+day[2:]+'_QC_met.dat'
                p_df=read_nsslmm(File,tstart,tend)
                masking = True #would you like a copy of this dataframe that is masked?
            except: 
                File= False #something has occured either the platform was not deployed or there was an error
                error_printing(e_test)
       
        elif p_type == 'UNL':
            try:
                File='/Volumes/Samsung_T5/Research/TORUS_Data/'+day+'/mesonets/UNL/UNL.'+pname+'.*'
                p_df = read_unlmm(File,tstart,tend)
                masking = False #would you like a copy of this dataframe that is masked?
            except: 
                File= False #something has occured either the platform was not deployed or there was an error
                error_printing(e_test)
        elif p_type == 'UAS':
            #To be filled in 
            File=False
            pass

        #If the data was successfully read in then add it to the dict
        if File != False:
            #pass the dataframe to a key item 'data' 
            pform[p_type][pname]['data']=p_df
            #assign a name for the data frame to its own dict key 
            pform[p_type][pname]['name']='{}_df'.format(pname)
            #attached said name to the pandas dataframe itself
            p_df.name='{}_df'.format(pname)
            
            if masking == True:
                #call the maskdata function 
                masked_df = maskdata(p_var, p_df)
                pform[p_type][pname]['masked_data']=masked_df
                pform[p_type][pname]['min']=masked_df.min()
                pform[p_type][pname]['max']=masked_df.max()
            elif masking == False: 
                #set max and min off of the unmasked dataframe
                pform[p_type][pname]['min']=p_df.min()[p_var]
                pform[p_type][pname]['max']=p_df.max()[p_var]
        testttt=pform[p_type][pname]
        print(pname,testttt['min'])

def find_by_key(data, target):
    for key, value in data.items():
        print(key)
        if isinstance(value, dict):
            for i in find_by_key(value, target):
                yield i 
        elif key == target:
            yield value 
for x in find_by_key(pform, "min"):
    print(x)

def gen_dict_extract(key, var):
    if hasattr(var,'iteritems'):
        for k, v in var.iteritems():
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result
# Convert Nested dictionary to Mapped Tuple 
# Using list comprehension + generator expression 
res = [(key, tuple(sub[key] for sub in pform.values())) for key in pform['NSSL']['FFld']] 
  
# printing result  
print("The grouped dictionary : " + str(res))  

#  # Convert Nested dictionary to Mapped Tuple
#  # Using defaultdict() + loop
#  res = defaultdict(tuple)
#  for key, val in test_dict.items():
#      for ele in val:
#          res[ele] += (val[ele], )
#
#  # printing result
#  print("The grouped dictionary : " + str(list(res.items()))
#
# Unnest single Key Nested Dictionary List 
# Using list comprehension 
#  res = {x : y[data_key] for idx in test_list for x, y in idx.items()}
testttt=pform['NSSL']['FFld']
print(testtttt['min'])

#  Unless I'm missing something you'd like to implement, I would go for a simple loop instead of using recursion:

def nested_get(input_dict, nested_key):
    internal_dict_value = input_dict
    for k in nested_key:
        internal_dict_value = internal_dict_value.get(k, None)
        if internal_dict_value is None:
            return None
    return internal_dict_value

print(nested_get({"a":{"b":{"c":1}}},["a","b","c"])) #1
print(nested_get({"a":{"bar":{"c":1}}},["a","b","c"])) #None


from jsonpath_rw import parse
def nested_get(d, path):
    result = parse(".".join(path)).find(d)
    return result[0].value if result else None


NSSLMM_df, UNLMM_df, max_array, min_array = [], [], [], []
for MM in ['FFld','LIDR','Prb1','Prb2','WinS','CoMeT1','CoMeT2','CoMeT3','UAS']:
   
    #load the NSSL mesonets 
    if (MM in ['FFld','LIDR','Prb1','Prb2','WinS']):
        try: 
            MMfile='/Users/severe2/Research/TORUS_data/'+day+'/mesonets/NSSL/'+MM+'_'+day[2:]+'_QC_met.dat'
        
            if MM == 'FFld':
                FFld_df = read_nsslmm(MMfile,tstart,tend)
                FFld_df.name='FFld'
                NSSLMM_df.append(FFld_df)
                masked_df = maskdata(p_var, FFld_df, mask=True)
                max_val, min_val= masked_df.max(), masked_df.min()
            elif MM == 'LIDR':
                LIDR_df = read_nsslmm(MMfile,tstart,tend)
                LIDR_df.name='LIDR'
                NSSLMM_df.append(LIDR_df)
                masked_df = maskdata(p_var, LIDR_df, mask=True)
                max_val, min_val= masked_df.max(), masked_df.min()
            elif MM == 'Prb1':
                Prb1_df = read_nsslmm(MMfile,tstart,tend)
                Prb1_df.name='Prb1'
                NSSLMM_df.append(Prb1_df)
                masked_df = maskdata(p_var, Prb1_df, mask=True)
                max_val, min_val= masked_df.max(), masked_df.min()
            elif MM == 'Prb2':
                Prb2_df = read_nsslmm(MMfile,tstart,tend)
                Prb2_df.name='Prb2'
                NSSLMM_df.append(Prb2_df)
                masked_df = maskdata(p_var, Prb2_df, mask=True)
                max_val, min_val= masked_df.max(), masked_df.min()
            elif MM == 'WinS':
                WinS_df = read_nsslmm(MMfile,tstart,tend)
                WinS_df.name='WinS'
                NSSLMM_df.append(WinS_df)
                masked_df = maskdata(p_var, WinS_df, mask=True)
                max_val, min_val= masked_df.max(), masked_df.min()
                #max_val, min_val= WinS_df[p_var].max(), WinS_df[p_var].min()
        except: error_printing(e_test)

    #load the UNL MM files 
    elif (MM in ['CoMeT1', 'CoMeT2', 'CoMeT3']):
        try: 
            MMfile='/Volumes/Samsung_T5/Research/TORUS_Data/'+day+'/mesonets/UNL/UNL.'+MM+'.*'
            
            if MM == 'CoMeT1':
                CoMeT1_df = read_unlmm(MMfile,tstart,tend)
                CoMeT1_df.name='CoMeT1'
                UNLMM_df.append(CoMeT1_df)
                masked_df = maskdata(p_var, CoMeT1_df, mask=True)
                max_val, min_val= masked_df.max(), masked_df.min()
            elif MM == 'CoMeT2':
                CoMeT2_df = read_unlmm(MMfile,tstart,tend)
                CoMeT2_df.name='CoMeT2'
                UNLMM_df.append(CoMeT2_df)
                masked_df = maskdata(p_var, CoMeT2_df, mask=True)
                max_val, min_val= masked_df.max(), masked_df.min()
            elif MM == 'CoMeT3':
                CoMeT3_df = read_unlmm(MMfile,tstart,tend)
                CoMeT3_df.name='CoMeT3'
                UNLMM_df.append(CoMeT3_df)
                masked_df = maskdata(p_var, CoMeT3_df, mask=True)
                max_val, min_val= masked_df.max(), masked_df.min()
        except: error_printing(e_test)

        print('********')
    #load the UAS files
    elif (MM in ['UAS','UAS_fillers']):
        UAS_df='Will come back and fill in'

    max_array.append(max_val)
    min_array.append(min_val)

#max an min values that will plot
#if p_var =='Thetae':
#    globalamin, globalamax=325, 360
#elif p_var =='Thetav':
#    globalamin, globalamax=295, 315

#Alternate global variable max/min for plotting purposes
globalamax = np.max(max_array)
globalamin = np.min(min_array)



#######################################################################################################
#################
# Create Plots ##
#################

#Create a Dictionary to contain the info relating to the plotting variable
var={
    'name': p_var,
    'min': globalamin,
    'max': globalamax
} 

#Create a Dictionary containing the insitu platform dataframes
pform={
    'NSSL': NSSLMM_df,
    'UNL': UNLMM_df,
    'UAS': UAS_df
}

if radar_plotting == True:
    print('\n Yes Print Radar \n')
    
    #get radar files
    ka1_files =sorted(glob.glob(filesys+'TORUS_Data/'+day+'/radar/TTUKa/netcdf/ka1/dealiased_*'))
    ka2_files =sorted(glob.glob(filesys+'TORUS_Data/'+day+'/radar/TTUKa/netcdf/ka2/dealiased_*'))
    fil= ka1_files + ka2_files
    
    ####dummy plot for scatter
    cmap=cmocean.cm.curl
    Z = [[0,0],[0,0]]
    levels = np.arange(globalamin,globalamax+1,1)
    CS3 = plt.contourf(Z, levels, cmap=cmap)#,transform=datacrs)
    plt.clf()

    
    #### what azimuth scan do you want to plot
    p_azimuth= 1.0

    #open radar data
    for thefile in fil[:]:
        print(str(thefile))
        radar = pyart.io.read(thefile)
        
        if radar.scan_type == 'ppi':
        
            n_swps = radar.nsweeps
            for swp_id in range(n_swps):
                plotter = pyart.graph.RadarDisplay(radar)
                azimuth = radar.fixed_angle['data'][swp_id]

                if np.around(azimuth, decimals=1) == p_azimuth:
                    print(str(thefile))
                    print("Producing Radar Plot")
                    #assign radar feilds and masking
                    det_radar_feilds(radar)

                    #get radar location data
                    currentscantime = datetime.strptime(radar.time['units'][14:-1], "%Y-%m-%dT%H:%M:%S")
                    try:
                        klat1, klon1, head1, rhib1, rhie1= det_radar_deps(currentscantime,'ka1',e_test)
                    except: error_printing(e_test)
                    try:
                        klat2, klon2, head2, rhib2, rhie2= det_radar_deps(currentscantime,'ka2',e_test)
                    except: error_printing(e_test)

                    #create dictionary containing the radar deployment data 
                    ka_info = {
                        'Main_Radar': 'Place Holder',
                        'Rfile': radar,
                        'Time': currentscantime,
                        'Swp_ID': swp_id,
                        'Azimuth': azimuth,
                        'KA1' : {
                            'lat': klat1,
                            'lon': klon1,
                            'head': head1,
                            'rhib': rhib1,
                            'rhie': rhie1,
                        },
                        'KA2' : {
                            'lat': klat2,
                            'lon': klon2,
                            'head': head2,
                            'rhib': rhib2,
                            'rhie': rhie2,
                        }
                    }

                    #Det which radar is the main plotting radar and redefine the 'main_radar' var for that radar's dict
                    if thefile.split("/")[-1][10:13] == 'Ka1':
                        ka_info['Main_Radar']= 'KA1'
                    elif thefile.split("/")[-1][10:13] == 'Ka2':
                        ka_info['Main_Radar']= 'KA2'
                    else:
                        print('What radar are we using?')
                    
                    #proceed to plot the radar
                    ppiplot(ka_info, var, pform, r_only, filesys, day, print_long, e_test)

                else:
                    print("Scan does not match the azimuth currently being plotted")
                    print("Currently plotting: azimuth= " + str(p_azimuth))
                    print("This scan: azimuth= "+ str(azimuth))
                    print('**********************************')
      
        elif radar.scan_type == 'rhi':
            print("This scan is an RHI")
            print('**********************************')
    
    print("all finished")

else:
#Only plotting timeseries (this code isn't fully fleshed out but in theroy this code is built in such a way to allow for this)
    print("no radar plotted")
    time_series(filesys,day,fig,p_var,radar_sub=False)
    #print(str(filesys+'TORUS_Data/'+day+'/mesonets/NSSL/*.nc'))

    fig.savefig('test2.png')
    plt.close()
