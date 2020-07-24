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
from read_platforms import read_nsslmm, read_unlmm
import json
import nexradaws
import pytz
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
conn = nexradaws.NexradAwsInterface()

##########################################################################
## VARIABLES
############
day='20190517' #'YYYYMMDD'
radar_plotting= True #would you like to plot radar images? (set to False to only plot timeseries)
r_only= False #Set to True for only radar as output, Set to False for Radar + Timeseries
p_var = "Thetav" #which var to plot (current options; Thetae, Thetav)
probe_of_interest='Prb1'
filesys='/Users/severe2/Research/'
temploc='/Volumes/Samsung_T5/Research/TORUS_Data/'

#Troubleshooting y/n
####################
#Set to True to print statements when entering/leaving definitions
print_long = False # True/False variable

#Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting 
  ##there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = False# True/False variable

def error_printing(e_test):
    ''' Basically I got sick of removing/replacing the try statements while troubleshooting 
    '''
    if e_test == True:
        e = sys.exc_info()
        traceback.print_tb(e[-1])
        #print( "<p>Error: %s</p>" % e )
        print("Error:", e[:-2])
        print(' ')
    else: pass
    return

##########################################################################
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

def get_WSR_from_AWS(start, end, radar_id,temploc):
    '''
    Retrieve the NEXRAD files that fall within a timerange for a specified radar site from the AWS server
    ----------
    INPUTS
    radar_id : string
        four letter radar designation
    start: datetime 
        start of the desired timerange
    end: datetime
        end of the desired timerange
    temploc: string
        location for the downloaded radarfiles
    -------
    RETURN
    radar_list : Py-ART Radar Objects
    '''
    #Determine the radar scans that fall within the time range for a given radar site
    scans = conn.get_avail_scans_in_range(start, end, radar_id)
    print("There are {} scans available between {} and {}\n".format(len(scans), start, end))
    
    #download the files that were identified
    #results = conn.download(scans[0:4], filesys+'TORUS_Data/'+day+'/radar/Nexrad/Nexrad_files/', keep_aws_folders=False)
    #results = conn.download(scans, filesys+'TORUS_Data/'+day+'/radar/Nexrad/Nexrad_files/', keep_aws_folders=False)
    results = conn.download(scans[0:4], temploc+day+'/radar/Nexrad/Nexrad_files/', keep_aws_folders=False)
    print("{} downloads failed: {}\n".format(results.failed_count,results.failed))
    #print("Results.iter_success : {}\n".format(results.iter_success()))
    
    #open the downloaded files as pyart objects
    radar_list=[]
    for i,scanfile in enumerate(results.iter_success()):
        #print("[{}] open_pyart, scan file_name = {}\n".format(i, scanfile.filename))
        try:
            radar_list.append(scanfile.open_pyart())
        except: 
            print("Failed to convert file: ",scanfile.filename)
    return radar_list

###########################
## PLOTTING DEFINTIONS  ###
###########################
def maskdata(p_var, platform_file, mask=True):
    platform_unmasked= platform_file[p_var].values
    
    if mask== False:
        platform_data= platform_unmasked
    elif mask== True:
        platform_name= str(platform_file.name)
        if (platform_name in ['FFld','WinS','LIDR','Prb1','Prb2']):
            platform_data= np.ma.masked_where(platform_file['qc_flag'].values>0, platform_unmasked)
        elif (platform_name in ['CoMeT1','CoMeT2','CoMeT3']):
            #for now
            platform_data= platform_unmasked
        elif (platform_name in ['Insert UAS filenames here']):
            print("will be filled in")
        else:
            print("What platform are you trying to use?")
    
    return platform_data

#* * * * * 
def getLocation(file, currentscantime, p_var, offsetkm, given_bearing= False):
    ''' 
    This definition has two functions:
        1) If no bearing is specified it will return the max and min lat/lons 
            to form a square surrounding the point indicated by lat1,lon1 by x km. 
            In this scenario end_lat and end_lon will be returned as nan's. 
        2) If a bearing is given then the defintion will return one set of lat/lon values 
             (end_lat and end_lon) which is the location x km away from the point lat1, lon1
             if an observer moved in a straight line the direction the bearing indicated. In 
             this scenario the four  min/max_lat/lon will be returned as nan's. 
    '''
    #determine the position of the probe at the time of radar scan
    p_sub, lon_sub, lat_sub, u_sub, v_sub, p_deploy = grab_platform_subset(file, currentscantime, p_var)
    current_lon=lon_sub[int(len(lon_sub)/2)]
    current_lat=lat_sub[int(len(lat_sub)/2)]

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

            if brng == 0:
                max_lat=new_lat
            elif brng == 90:
                max_lon= new_lon
            elif brng == 180: 
                min_lat= new_lat
            elif brng == 270: 
                min_lon=new_lon 
        end_lat, end_lon = np.nan, np.nan 

    else: 
        #if a bearing is provided 
        bearing = (given_bearing/ 90.)* np.pi / 2.

        new_lat = np.arcsin(np.sin(lat1) * np.cos(offsetkm/R) + np.cos(lat1) * np.sin(offsetkm/R) * np.cos(bearing))
        new_lon = lon1 + np.arctan2(np.sin(bearing)*np.sin(offsetkm/R)*np.cos(lat1),np.cos(offsetkm/R)-np.sin(lat1)*np.sin(new_lat))
        end_lon = 180.0 * new_lon / np.pi
        end_lat = 180.0 * new_lat / np.pi
        min_lat, min_lon, max_lat, max_lon = np.nan, np.nan, np.nan, np.nan   

    return min_lat, max_lat, min_lon, max_lon, end_lat, end_lon 

#* * * * * 
def legend_maker(p_name,m_style,m_color,leg_str,legend_elements):
    if print_long== True: print('made it into legend_maker')
    
    if (p_name in ['Ka1','Ka2']): #the legend entries for the KA radars
        legend_entry=Line2D([], [], marker=m_style, markeredgecolor='black',markeredgewidth=3,label=leg_str,markerfacecolor=m_color, markersize=26)
    else:
        legend_entry=Line2D([], [], marker=m_style, markeredgecolor=m_color,markeredgewidth=3,label=leg_str, markersize=26,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')])
    
    legend_elements=np.append(legend_elements,legend_entry)
    
    if print_long== True: print('made it through legend maker')
    return legend_elements

#* * * * * 
def platform_attr(p_file,legend_elements,radar_m=False,r_s=False):
    if print_long== True: print('Made it into platform_attr')

    #assign the atributes for non-radar markers
    if radar_m == False:
        if p_file.name == "Prb1":
            m_style, m_color, l_color, leg_str= '1','xkcd:lightblue','steelblue','Prb1'
        elif p_file.name == "Prb2":
            m_style, m_color, l_color, leg_str= '1','xkcd:watermelon','xkcd:dusty red','Prb2'
        elif p_file.name == "FFld":
            m_style, m_color, l_color, leg_str= '1','xkcd:bubblegum pink','xkcd:pig pink','FFld'
        elif p_file.name == "LIDR":
            m_style, m_color, l_color, leg_str= '1','xkcd:pastel purple','xkcd:light plum','LIDR'
        elif p_file.name == "WinS":
            m_style, m_color, l_color, leg_str= '1','xkcd:peach','xkcd:dark peach','WinS'
        elif p_file.name == "CoMeT1":
            m_style, m_color, l_color, leg_str= '1','brown','brown','CoMeT1'
        elif p_file.name == "CoMeT2":
            m_style, m_color, l_color, leg_str= '1','yellow','yellow','CoMeT2'
        elif p_file.name == "CoMeT3":
            m_style, m_color, l_color, leg_str= '1','black','black','CoMeT3'
        p_name = p_file.name

    #assign the atributes for the radar markers
    elif radar_m == True:
        if p_file == 'Ka2':
            m_style, m_color, l_color, leg_str= '8','xkcd:crimson','xkcd:crimson','Ka2'
        elif p_file == 'Ka1':
            m_style, m_color, l_color, leg_str= '8','mediumseagreen','mediumseagreen','Ka1'
        p_name = p_file

    if r_s == True: #only add to the legend if the platform_attr def was called while making the radar subplots
        legend_elements=legend_maker(p_name,m_style,m_color,leg_str,legend_elements)
    
    if print_long== True: print('Made it through platform_attr')
    return m_style, m_color, l_color, leg_str, legend_elements

# * * * * * *
def plot_colourline(x,y,c,cmap,ax,datacrs,color='k',amin=None,amax=None):
    if amin:
        pass
    else:
        amin=np.nanmin(c)
    if amax:
        pass
    else:
        amax=np.nanmax(c)
        
    c = cmap((c-amin)/(amax-amin))
    ax = plt.gca()
    for i in np.arange(len(x)-1): #This is the border of the colorline
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=color, linewidth=10.5, transform=datacrs)
    for i in np.arange(len(x)-1): #This is the colorramp colorline
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=c[i], linewidth=7.5, transform=datacrs)
    
    return

# * * * * * *
def grab_platform_subset(p_df, scan_time, p_var):
    ''' Grabs the observed thermo data +/-5 minutes around radar scan time
        (does not automatically control the gray vertical box on timeseries)
    '''
    aaa = p_df.loc[(p_df['datetime']>=scan_time-dt.timedelta(seconds=300))]
    df_sub = aaa.loc[(aaa['datetime']<=scan_time+dt.timedelta(seconds=300))]

    p_sub = df_sub[p_var].values
    U_sub = df_sub['U'].values
    V_sub = df_sub['V'].values
    lat_sub= df_sub['lat'].values
    lon_sub= df_sub['lon'].values
    
    # test to ensure that there is valid data in the subrange (aka the platform was deployed during the time of radarscan)
    try:
        # p_test=df_sub.loc[df_sub['datetime']==scan_time
        p_test= df_sub.iloc[1]
        p_deploy = True
    except:
        print('The platform was not deployed at this time')
        p_deploy = False
        error_printing(e_test)

    return p_sub, lon_sub, lat_sub, U_sub, V_sub, p_deploy

# * * * * * *
def platform_plot(file,radartime,ax_n,color,m_color,p_var,e_test,labelbias=(0,0)):
    if print_long==True: print('made it into platform_plot')

    #grab the subset of data of +- interval around radar scan
    p_sub, lon_sub, lat_sub, u_sub, v_sub, p_deploy = grab_platform_subset(file, radartime, p_var)
    
    if p_deploy == True:
        #plot data for reflectivity plot
        plot_colourline(lon_sub,lat_sub,p_sub,cmocean.cm.curl,ax_n,ccrs.PlateCarree(),color=color,amin=globalamin,amax=globalamax)
  
        #d2=plt.text(lon_sub[int(len(lon_sub)/2)]+labelbias[0], lat_sub[int(len(lat_sub)/2)]+labelbias[1], file[58:62],#transform=ccrs.PlateCarree(), fontsize=20, zorder=9, path_effects=[patheffects.withstroke(linewidth=4,foreground=color)])
        ax_n.plot(lon_sub[int(len(lon_sub)/2)], lat_sub[int(len(lat_sub)/2)], transform=ccrs.PlateCarree(), marker='1',markersize=20, markeredgewidth='3',color=m_color,path_effects=[patheffects.withStroke(linewidth=12,foreground='k')], zorder=10)
  
        ax_n.plot(lon_sub[int(len(lon_sub))-1], lat_sub[int(len(lat_sub))-1], transform=ccrs.PlateCarree(), marker='.',markersize=10, markeredgewidth='3',color='k',zorder=9) 
  
        #plot windbarbs
        try:
            stationplot = metpy.plots.stationplot(ax_n, lon_sub[::30], lat_sub[::30], clip_on=true, transform=ccrs.PlateCarree())
            stationplot.plot_barb(u_sub[::30], v_sub[::30],length=7)
        except: error_printing(e_test)

    if print_long==True: print('made it through platform_plot')
    return



##  #  #  #  #  #  #  #  #  #  #  #  #  # # # #  #  #  #  #  #  #  #  #  #  # #  #  #  # # #  #  #  #  #  #  #  #  # #  #
##  #  #  #  #  #  #  #  #  #  #  #  #  # # # #  #  #  #  #  #  #  #  #  #  # #  #  #  # # #  #  #  #  #  #  #  #  # #  #
def ppiplot(r_only,radar_list,radar,filesys,day,globalamin, globalamax,p_var,p_of_int,CS3,e_test):
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
    post=radar_subplots('refl',swp,fig,display,currentscantime,globalamin,globalamax,p_var,p_of_int,gs,e_test)
    post=radar_subplots('vel',swp,fig,display,currentscantime,globalamin,globalamax,p_var,p_of_int,gs,e_test)
    if r_only == False:
        time_series(filesys,day,fig,currentscantime,globalamin, globalamax, p_var,gs,radar_sub=True)
  
    ## Plot platform colorbar
    if p_var == "Thetae":
        c_lab= "Equivalent Potential Temp [K]"
    elif p_var == "Thetav":
        c_lab= "Virtual Potential Temp [K]"
    cbar_ax = plt.axes([.514, post.y0,.014, post.y1-post.y0])#left, bottom, width, height
    cbar=plt.colorbar(CS3,cax=cbar_ax, orientation='vertical', label=c_lab, ticks=MaxNLocator(integer=True))#,ticks=np.arange(globalamin,globalamax+1,2))

    ## Plot title
    title=r_site+' '+str(elev)+r'$^{\circ}$ PPI '+ fancy_date_string_utc
    plt.suptitle(title,y=.92)

    ## Finish plot 
    plt.savefig('/Users/severe2/Research/TORUS_Data/'+day+'/radar/Nexrad/plots/'+currentscantime.strftime('%m%d%H%M')+'_'+r_site+'_'+p_var+'.png' ,bbox_inches='tight',pad_inches=.3)
    plt.close()

    if print_long== True: print('~~~~~~~~~~~made it through ppiplot~~~~~~~~~~~~~~~~~~')
    print('Done Plotting')
    print('******************************************************************************************************')

    ## Makes a ding noise 
    print('\a')
    return

def radar_subplots(mom,swp,fig,display,currentscantime,globalamin,globalamax,p_var,p_of_int,gs,e_test):
    if print_long== True: print('~~~~~~~~~~~made it into radar_subplots~~~~~~~~~~~~~~')

    ## SET UP PLOTTING CONTROLS
    NSSLmm, NEBmm, UASd = True, True, False #which other platforms do you wish to plot 
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
    ymin,ymax,xmin,xmax,a,b=getLocation(p_of_int,currentscantime,p_var,offsetkm=40)    
   
    # SET UP SUBPLOTS
    ax_n = fig.add_subplot(gs[row,col], projection=display.grid_projection)
    ax_n.text(.5,-.065,p_title, transform = ax_n.transAxes, horizontalalignment = 'center', fontsize = 40) #the radar subplot titles
    #  display.plot_ppi_map(feild, swp_id, title_flag=False,colorbar_flag=False,cmap=c_scale,ax=ax_n,vmin=vminb,vmax=vmaxb,embelish=False)
    display.plot_ppi_map(feild, swp_id, title_flag=False, colorbar_flag=False, cmap=c_scale, ax=ax_n,vmin=vminb,vmax=vmaxb,min_lon=xmin,max_lon=xmax,min_lat=ymin,max_lat=ymax,embelish=False)
    
    ## PLOT OTHER PLATFORMS             
    if NSSLmm == True:
        for NSSLMM in NSSLMM_df:
            if print_long== True: print(NSSLMM.name)
            m_style,m_color,l_color,leg_str,legend_elements=platform_attr(NSSLMM,legend_elements,r_s=True)
            platform_plot(NSSLMM,currentscantime,ax_n,'xkcd:light grey',m_color,p_var,e_test, labelbias=(0,0))
    if NEBmm == True:
        for UNLMM in UNLMM_df:
            if print_long== True: print(UNLMM.name)
            m_style,m_color,l_color,leg_str,legend_elements=platform_attr(UNLMM,legend_elements,r_s=True)
            platform_plot(UNLMM,currentscantime,ax_n,'xkcd:light grey',m_color,p_var,e_test, labelbias=(0,0))
    if UASd == True:
        for UAS in UAS_files:
            m_style,m_color,l_color,leg_str,legend_elements=platform_attr(UAS,legend_elements,r_s=True)
            platform_plot(UAS,currentscantime,ax_n,'xkcd:very pale green',p_var,labelbias=(0,0.01))
 
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

def time_series(filesys,day,fig,currentscantime,globalamin,globalamax,p_var,gs,radar_sub=False):
    if print_long== True: print('~~~~~~~~~~~Made it into time_series~~~~~~~~~~~~~~~~~')
     
    if radar_sub == True:
        vline=currentscantime 
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
                if print_long== True: print(str(platform_file.name))
                
                ## Set up line colors and legend labels 
                b_array=[] #blank array
                m_style,m_color,l_color,leg_str,bb_array = platform_attr(platform_file,b_array) #get the attributes (color, label, etc)

                ## Should ploting data be masked?
                # if you only want to mask certain platforms set up an if-statement to assign the "masking" variable
                masking = True #True/False variable
                plot_data= maskdata(p_var,platform_file,masking)

                ## Plot
                ax_n.plot(platform_file['datetime'],plot_data,linewidth=3,color=l_color,label=leg_str) #assigning label= is what allows the legend to work  
                    
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
        ax_n.vlines(vline,globalamin-1,globalamax+1,colors='r',linewidth=4, alpha=.5)
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
            
            if MM == probe_of_interest:
                p_of_int=NSSLMM_df[-1]
        
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

            if MM == probe_of_interest:
                p_of_int=UNLMM_df[-1]

        except: error_printing(e_test)

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
#Dummy plot for scatter
# Make a scale covering max and min
cmap=cmocean.cm.curl
Z = [[0,0],[0,0]]
levels = np.arange(globalamin,globalamax+1,1)
CS3 = plt.contourf(Z, levels, cmap=cmap)#,transform=datacrs)
plt.clf() # clear current figure
# ************************************

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
    radar_list=get_WSR_from_AWS(ts[r],te[r],r,temploc)
    for radar in radar_list:
        #print(radar.info())
        #testing_plots(radar,i,j,r_only,globalamin,globalamax,p_var,e_test)
        #plot_with_WSR(radar,i,j)

        #why am i calling radar and radarlist to this function when I am in a loop that is interating through these to values?
        ppiplot(r_only, radar_list,radar, filesys, day, globalamin,globalamax, p_var, p_of_int,CS3, e_test)
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
