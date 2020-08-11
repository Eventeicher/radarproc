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
print_long=config.print_long
e_test=config.e_test
p_var=config.p_var
filesys=config.filesys
temploc=config.temploc
day=config.day

## Read in defns that I have stored in another file (for ease of use/consistancy accross multiple scripts)
from shared_defns import Add_to_DATA, pform_names, error_printing

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
    ''' This definition has two functions:
            1) If no bearing is specified it will return a namedtuple containing the max/min lat/lons
                to form a square surrounding the point indicated by lat1,lon1 by x km.
            2) If a bearing is given then the defintion will return one set of lat/lon values
                 (end_lat and end_lon) which is the location x km away from the point lat1, lon1
                 if an observer moved in a straight line the direction the bearing indicated.
    ----
    INPUTS: current_lon & current_lat: the starting position
            offsetkm: distance traveled from the starting point to the new locations
            given_bearing: True/False, are you given a specified direction of "travel"
    '''
    lat1 = current_lat * np.pi / 180.0
    lon1 = current_lon * np.pi / 180.0
    R = 6378.1 #earth radius (R = ~ 3959 MilesR = 3959)

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
    latArray,lonArray = [], []
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
                ang= ka_info[ka]['head'] + j
                if ang > 360.: ang=int(ang-360.)

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
        INPUTS: x & y: the loction of the probe at a given time
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
        except: error_printing(e_test)

    elif p_deploy == False:
        if print_long==True: print('The platform was not deployed at this time')

    if print_long== True: print('made it through platform_plot')
    return



##  #  #  #  #  #  #  #  #  #  #  #  #  # # # #  #  #  #  #  #  #  #  #  #  # #  #  #  # # #  #  #  #  #  #  #  #  # #  #
##  #  #  #  #  #  #  #  #  #  #  #  #  # # # #  #  #  #  #  #  #  #  #  #  # #  #  #  # # #  #  #  #  #  #  #  #  # #  #
def ppiplot(ka_info, var, pform, filesys, day, print_long, e_test):
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
    if t_plotting == False:
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
    if t_plotting == True: time_series(fig, gs, ka_info, var, pform, print_long, e_test)

    ## Plot platform colorbar
    if var['name'] == "Thetae": c_lab= "Equivalent Potential Temp [K]"
    elif var['name'] == "Thetav": c_lab= "Virtual Potential Temp [K]"
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
    print('Done Plotting \n ***************************************************************************************************')

    ## Makes a ding noise
    print('\a')
    return

# * * * * * *  *
def radar_subplots(mom, fig, gs, display, ka_info, var, pform, filesys, print_long, e_test):
    ''' Plots each of the radar subplots and calls for the markers to be plotted for all the additional platforms
    ----
    INPUTS: mom: str, indicates which radar moment you would like to plot (ie Reflectivity, Velocity, Spectrum Width etc)
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
    if rhi_ring == True: rhi_spokes_rings(ka_info, display, print_long, e_test) #plot the spokes for radars

    ## PLOT RADAR MARKERS
    for ka in ['KA1', 'KA2']:
        try:
            if np.logical_and(ka_info[ka]['lat']>box.ymin,np.logical_and(ka_info[ka]['lat']<box.ymax,np.logical_and(ka_info[ka]['lon']>box.xmin,ka_info[ka]['lon']<box.xmax))):
                P_Attr,legend_elements=pform_attr(ka,print_long,r_s=True,radar_m=True,l_array=legend_elements)

                ax_n.plot(ka_info[ka]['lon'],ka_info[ka]['lat'],marker=P_Attr['m_style'], transform=ccrs.PlateCarree(),color=P_Attr['m_color'],markersize=18,markeredgewidth=5,path_effects=[PathEffects.withStroke(linewidth=15,foreground='k')],
                        zorder=10)
                #d2=ax_n.text(ka_info[ka]['lon']+0.006, ka_info[ka]['lat']-0.011,'Ka2',transform=ccrs.PlateCarree(), zorder=10, path_effects=[PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')]) #textlabel on plot
        except: error_printing(e_test)

    ## PLOT INSITU TORUS PLATFORMS
    if NSSLmm == True:
        for NSSLMM in NSSLMM_df:
            if print_long== True: print(NSSLMM.name)
            P_Attr,legend_elements=pform_attr(NSSLMM, print_long, r_s=True, l_array=legend_elements)
            platform_plot(NSSLMM, ka_info['Time'], var, P_Attr, ax_n, print_long, e_test)
    if NEBmm == True:
        for UNLMM in UNLMM_df:
            if print_long== True: print(UNLMM.name)
            P_Attr,legend_elements=pform_attr(UNLMM, print_long, r_s=True, l_array=legend_elements)
            platform_plot(UNLMM, ka_info['Time'], var, P_Attr, ax_n, print_long, e_test)
    if UASd == True:
        for UAS in UAS_files:
            P_Attr,legend_elements=pform_attr(UAS, print_long, r_s=True, l_array=legend_elements)
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
    else: pass #this means you are currently making the right subplot

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
def time_series(fig, gs, ka_info, var, pform, print_long, e_test, r_plotting=True):
    ''' Plot the time series of p_var from the various instrument platforms
    ----
    INPUTS; fig & gs: relating to the plot layout
            ka_info, var, & pform: dictionarys, as described in the ppiplot comment
                **** NOTE: I need to come back and do more with pform****
            r_plotting: True/False, If False will only plot timeseries (no radar).. In theroy have not actually done this yet
    '''
    if print_long== True: print('~~~~~~~~~~~Made it into time_series~~~~~~~~~~~~~~~~~')

    if r_plotting == True: ax_n=fig.add_subplot(gs[1,:])
    else: fig, ax_n= plt.subplots() #will only plot the time series (not the radar subplots)

    ## MAKE THE TIMESERIES
    for platform_df in [NSSLMM_df,UNLMM_df,UAS_df]:
        #can remove this if/else statement once you fill in the UAS portion
        if platform_df == UAS_df:
            # Need to come back and fill in UAS
            pass
        else:
            for platform_file in platform_df:
                if print_long== True: print(str(platform_file.name))
                ## Set up line colors and legend labels
                P_Attr = pform_attr(platform_file, print_long)

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
    if var['name'] == "Thetae": ylab = "Equivalent Potential Temp [K]"
    elif var['name'] == "Thetav": ylab = "Virtual Potential Temp [K]"
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
    for line in leg.get_lines(): line.set_linewidth(12)

    if print_long== True: print('~~~~~~~~~~~Made it through time_series~~~~~~~~~~~~~~')
    return


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
    ka1_files =sorted(glob.glob(filesys+'TORUS_Data/'+day+'/radar/TTUKa/netcdf/ka1/dealiased_*'))
    ka2_files =sorted(glob.glob(filesys+'TORUS_Data/'+day+'/radar/TTUKa/netcdf/ka2/dealiased_*'))
    fil= ka1_files + ka2_files

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
                    #  Establish locations for other radars (if deployed)
                    Data, subset_pnames = Add_to_DATA('RADAR',Data,subset_pnames,print_long)
       
                    ## Proceed to plot the radar
                    ##### + + + + + + + + + + + +
                    #  ppiplot(Data, print_long, e_test)

                ## If the scans azimuth does not match the one we are plotting
                else: print("Scan does not match the azimuth currently being plotted \n Currently plotting: azimuth= " 
                             + str(config.p_azimuth) +"\n This scan: azimuth= "+ str(azimuth)+'\n**********************************')
print(Data)
print(subset_pnames)
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
