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
import matplotlib as mpl
import glob
import os
import matplotlib
import matplotlib.cm as cm
import matplotlib.patheffects as pe
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
from read_nssl import read_nsslmm
import json
import nexradaws
import pytz
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
conn = nexradaws.NexradAwsInterface()


## VARIABLES
############
day='20190517' #'YYYYMMDD'
filesys='/Users/severe2/Research/'
p_var= 'Thetav'
probe_of_interest='Prb1'
temploc='/Users/severe2/Research/TORUS_Data/'
r_only= False #Set to True for only radar as output, Set to False for Radar + Timeseries

#crop the start time to a time you specify (comment out to do the whole day) 
#could also set up a tend var in same format
tstart = dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),15,0,0)
tend = None

#Troubleshooting y/n
####################
#Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting 
  ##there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = True #True/False variable

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

## Functions to get radar data
################################
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
    results = conn.download(scans, filesys+'TORUS_Data/'+day+'/radar/Nexrad/Nexrad_files/', keep_aws_folders=False)
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

# * * * * * *
def maskdata(p_var, NSSLMM):
    p_um=NSSLMM[p_var].values
    p_m=np.ma.masked_where(NSSLMM['qc_flag'].values>0, p_um)
    return p_m, p_um

def getLocation(file, currentscantime,p_var,offset=0):
    ''' Determine the location surrounding a given probe for a given offset
    '''
    #determine the position of the probe at the time of radar scan
    p_sub, lon_sub, lat_sub, u_sub, v_sub, p_deploy = grab_platform_subset(file, currentscantime, p_var)
    current_lon=lon_sub[int(len(lon_sub)/2)]
    current_lat=lat_sub[int(len(lat_sub)/2)]

    lat_min,lat_max = current_lat-offset, current_lat+offset
    lon_min,lon_max = current_lon-offset, current_lon+offset
    
    return lat_min, lat_max, lon_min, lon_max

def legend_maker(p_name,m_style,m_color,leg_str,legend_elements):
    print('made it into legend_maker')
    if (p_name in ['ka1','ka2']): #the legend entries for the ka radars
        legend_entry=Line2D([], [], marker=m_style, markeredgecolor='black',markeredgewidth=3,label=leg_str,markerfacecolor=m_color, markersize=26)
    else:
        legend_entry=Line2D([], [], marker=m_style, markeredgecolor=m_color,markeredgewidth=3,label=leg_str, markersize=26,path_effects=[patheffects.withStroke(linewidth=12,foreground='k')])
    
    #  if m_style == '8': #the legend entries for the KA radars
        #  legend_entry=line2D([], [], marker=m_style, markeredgecolor='black',markeredgewidth=3,label=leg_str,markerfacecolor=m_color, markersize=26)
    #  elif m_style == '1':
        #  legend_entry=line2D([], [], marker=m_style, markeredgecolor=m_color,markeredgewidth=3,label=leg_str, markersize=26,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')])
    
    legend_elements=np.append(legend_elements,legend_entry)
    print('made it through legend maker')
    
    return legend_elements

def platform_attr(p_file,legend_elements,radar_m=False,r_s=False):
    print('Made it into platform_attr')
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
        p_name=p_file.name

    elif radar_m == True:
        if p_file == 'Ka2':
            m_style, m_color, l_color, leg_str= '8','xkcd:crimson','xkcd:crimson','Ka2'
        elif p_file == 'Ka1':
            m_style, m_color, l_color, leg_str= '8','mediumseagreen','mediumseagreen','Ka1'
        p_name=p_file

    if r_s == True: #only add to the legend if the platform_attr def was called while making the radar subplots
        legend_elements=legend_maker(p_name, m_style,m_color,leg_str,legend_elements)
    print('Made it through platform_attr')

    return m_style, m_color, l_color, leg_str, legend_elements

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

def platform_plot(file,radartime,ax_n,color,m_color,p_var,e_test,labelbias=(0,0)):
    print('made it into platform_plot')

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

    print('made it through platform_plot')
    
    return

# **********************************
def ppiplot(r_only,radar_list,radar,filesys,day,globalamin, globalamax,p_var,p_of_int,e_test,i,j):
    print('made it into ppiplot')
    
    SMALL_SIZE, MS_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 25, 30, 35, 50
    plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=MS_SIZE)       # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    sweep=0
    # Find the beginning location of designated sweep in data array
    index_at_start = radar.sweep_start_ray_index['data'][sweep]
    # Time at which scanning began
    time_at_start_of_radar = num2date(radar.time['data'][index_at_start],radar.time['units'])
    print(time_at_start_of_radar)
    # Convert time into fancy date string to use in title
    fancy_date_string_utc = time_at_start_of_radar.strftime('%Y-%m-%d %H:%M UTC')
    print(fancy_date_string_utc)
    print('Gathered Radar Data')
    
    ## Set up plot size
    if r_only == True:
        fig = plt.figure(figsize=(40,150),facecolor='white')
    else:
        fig = plt.figure(figsize=(30,20),facecolor='white')
    display = pyart.graph.RadarMapDisplay(radar)

    ## Bounding box, x km away from the radar in all directions (set number back to 21)
    #xmin, xmax = getLocation(klat,klon,270,21)[1], getLocation(klat,klon,90,21)[1]
    #ymin, ymax = getLocation(klat,klon,180,21)[0], getLocation(klat,klon,0,21)[0]

    ## Plot instuments on subplots
    radar_subplots('refl',fig,display,time_at_start_of_radar,globalamin,globalamax,p_var,p_of_int,e_test)
    radar_subplots('vel',fig,display,time_at_start_of_radar,globalamin,globalamax,p_var,p_of_int,e_test)
    if r_only == False:
        time_series(filesys,day,fig,time_at_start_of_radar,globalamin, globalamax, p_var,radar_sub=True)
    
    # Dummy plot for scatter
    # Make a scale covering max and min
    #  cmap=cmocean.cm.curl
    #  Z = [[0,0],[0,0]]
    #  levels = np.arange(globalamin,globalamax+1,1)
    #  print(Z)
    #  print(levels)
    #  print(cmap)
    #  CS3 = plt.contourf(Z, levels, cmap=cmap)#,transform=datacrs)
    #  plt.clf() # clear current figure
#
    # Plot colorbars
    #  if p_var == "Thetae":
        #  c_lab= "Equivalent Potential Temp [K]"
    #  elif p_var == "Thetav":
        #  c_lab= "Virtual Potential Temp [K]"
    #  cbar_ax = plt.axes([.5265,.5356,.014, 0.405])#left, bottom, width, height
    #  cbar=plt.colorbar(CS3,cax=cbar_ax, orientation='vertical', label=c_lab, ticks=MaxNLocator(integer=True))#,ticks=np.arange(globalamin,globalamax+1,2))
#
    ## Plot title
    title=fancy_date_string_utc
    plt.suptitle(title)

    ## Finish plot 
    plt.tight_layout()
    plt.savefig('/Users/severe2/Research/TORUS_Data/'+day+'/radar/Nexrad/plots/goal'+str(i)+'_'+str(j)+'.png')
    plt.close()

    print('made it through ppiplot')
    print('***********************************************************************************************************')
    
    return

def radar_subplots(mom,fig,display,currentscantime,globalamin,globalamax,p_var,p_of_int,e_test):
    print("made it into radar sub")
    ## SET UP PLOTTING CONTROLS
    NSSLmm, NEBmm, UASd = True, False, False #which other platforms do you wish to plot 
    country_roads, hwys, county_lines, state_lines = False, False, False, False #background plot features
    legend_elements=[]
    ###########################

    ## SET UP VARS FOR EACH RADAR MOMENTS
    if mom == 'refl':
        pos, swp_id, feild = 221, 0, 'reflectivity'
        c_scale, c_label = 'pyart_HomeyerRainbow', 'Radar Reflectivity [dbz]'
        #vminb, vmaxb = -30., 30.
        vminb, vmaxb = -10., 75.
        p_title, leg = 'Reflectivity', True
    elif mom == 'vel':
        pos, swp_id, feild = 222, 1, 'velocity'
        c_scale, c_label = 'pyart_balance', 'Velocity [m/s]'
        vminb, vmaxb = -40., 40.
        p_title, leg = 'Radial Velocity', False
    else:
        print("hey what just happened!\n")
        exit

    ## Bounding box, x km away from the probe of interest in all directions 
    ymin,ymax,xmin,xmax=getLocation(p_of_int,currentscantime,p_var,offset=.75)    
    
    # SET UP SUBPLOTS
    ax_n = fig.add_subplot(pos, projection=display.grid_projection)
    ax_n.text(.5,-.065,p_title, transform = ax_n.transAxes, horizontalalignment = 'center', fontsize = 40) #the radar subplot titles
    display.plot_ppi_map(feild, swp_id, title_flag=False, cmap=c_scale, ax=ax_n,vmin=vminb,vmax=vmaxb,min_lon=xmin,max_lon=xmax,min_lat=ymin,max_lat=ymax,colorbar_label= c_label,embelish=False) 
    
    ## PLOT OTHER PLATFORMS             
    if NSSLmm == True:
        for NSSLMM in NSSLMM_df:
            print(NSSLMM.name)
            m_style,m_color,l_color,leg_str,legend_elements=platform_attr(NSSLMM,legend_elements,r_s=True)
            platform_plot(NSSLMM,currentscantime,ax_n,'xkcd:light grey',m_color,p_var,e_test, labelbias=(0,0))
    if NEBmm == True:
        print('To be filled in')
    if UASd == True:
        for UAS in UAS_files:
            m_style,m_color,l_color,leg_str,legend_elements=platform_attr(UAS,legend_elements,r_s=True)
            platform_plot(UAS,currentscantime,ax_n,'xkcd:very pale green',p_var,labelbias=(0,0.01))
 
    
    ## SET UP LEGENDS
    l_list=legend_elements.tolist()
    if leg == True: #add legend for platform markers
        l=ax_n.legend(handles=l_list,loc='lower left', bbox_transform=ax_n.transAxes, bbox_to_anchor=(-0.45,0), handlelength=.1,title="Platforms",shadow=True,fancybox=True,ncol=1,edgecolor='black')
        l.get_title().set_fontweight('bold')
    else: #Set up an invisible legend in a jankey method to force the plots to be where I want them (#goodenoughforgovwork)
        ax_n.legend([],[],loc='lower left', bbox_transform=ax_n.transAxes, bbox_to_anchor=(1.3,1.017), handlelength=.25,frameon=False)
    print('made it through radar sub')

    return

def time_series(filesys,day,fig,currentscantime,globalamin,globalamax,p_var,radar_sub=False):
    print('Made it into time_series')
    if radar_sub == True:
        ax_n, var1='ax3', 212
        vline=currentscantime 
        ax_n=fig.add_subplot(var1)
    
    else: #will only plot the time series (not the radar subplots)
        ax_n='ax'
        fig, ax_n= plt.subplots()
    
    #NSSLMM_files= sorted(glob.glob(filesys+'TORUS_Data/'+day+'/mesonets/NSSL/*.nc'))
    for NSSLMM in NSSLMM_df:
        print(str(NSSLMM.name))
        b_array=[] #blank array
        m_style,m_color,l_color,leg_str,bb_array = platform_attr(NSSLMM,b_array) #get the attributes (color, label, etc)

        ## OPEN the file and pull out data we want
        if p_var == "Thetae":
            ylab = "Equivalent Potential Temp [K]"
        elif p_var == "Thetav":
            ylab = "Virtual Potential Temp [K]"

        ## Plot
        #should ploting data be masked?
        mask = True #True/False variable
        if mask == True:
            p_mask, p_unmasked = maskdata(p_var, NSSLMM)
            ax_n.plot(NSSLMM['datetime'],p_mask,linewidth=3,color=l_color,label=leg_str) #assigning label= is what allows the legend to work  
        elif mask == False:
            p_mask, p_unmasked = maskdata(p_var, NSSLMM)
            ax_n.plot(NSSLMM['datetime'],p_unmasked,linewidth=3,color=l_color,label=leg_str) #assigning label= is what allows the legend to work  
        
        ## Set up XY axes tick locations
        ax_n.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M')) #strip the date from the datetime object
        ax_n.xaxis.set_major_locator(mdates.AutoDateLocator()) 
        ax_n.xaxis.set_minor_locator(AutoMinorLocator(6)) # set up minor ticks (should be a multiple of ten intervals (ie 10,20,30... min spans)
        ax_n.yaxis.set_major_locator(MultipleLocator(5)) # set up major tick marks (this is set up to go by 5's will want to change for diff vars)
        ax_n.yaxis.set_minor_locator(AutoMinorLocator(5)) # set up minor ticks (this have it increment by 1's will want to change for diff vars)
    
    ## Set up axes formats (lables/ ticks etc)
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
    
    print('Made it through time_series')
    return

#*****************************
#****************************************

def plot_with_WSR(radar,i,j):
    # Loop through each radar volume to plot
    sweep=0
    # Find the beginning location of designated sweep in data array
    index_at_start = radar.sweep_start_ray_index['data'][sweep]
    # Time at which scanning began
    time_at_start_of_radar = num2date(radar.time['data'][index_at_start],
                                      radar.time['units'])
    fancy_date_string_utc = time_at_start_of_radar.strftime('%Y-%m-%d %H:%M UTC')

    # Grab the location data for both radars
    # ka1_loc = k_loc[k_loc['radar']=='ka1'].reset_index()
    # ka2_loc = k_loc[k_loc['radar']=='ka2'].reset_index()

    # Begin figure
    fig = plt.figure(figsize = [20,20])
    display = pyart.graph.RadarMapDisplay(radar)
    proj = display.grid_projection
    lat_0 = display.loc[0] # Use radar data extent to define plot region
    lon_0 = display.loc[1]
    lat_1 = 33.96

    # Plot reflectivity on the left
    ax1 = fig.add_subplot(1,2,1,projection=proj)
    plt.tight_layout()
    field='reflectivity'
    title = field+' \n'+fancy_date_string_utc#+fancy_date_string+' ('+fancy_date_string_utc+') '
    offset = 0.75
    #ax1.plot(klon1,klat1,color='r',markersize=200,transform=proj)
    ax1.plot(transform=proj)
    display.plot_ppi_map(field, 0, colorbar_flag=True,title=title, ax=ax1,vmin=-10, vmax=75)#min_lon=klon1-offset,max_lon=klon1+offset,
            #min_lat=klat1-offset,max_lat=klat1+offset)
    
    ax2= fig.add_subplot(1,2,2,projection=proj)
    feild='velocity'
    title = field+' \n'+fancy_date_string_utc#+fancy_date_string+' ('+fancy_date_string_utc+') '
    offset = 0.75
    norm, cmap = ctables.registry.get_with_steps('NWSVelocity', 16, 16)
    display.plot_ppi_map(field, 0, colorbar_flag=True,title=title, ax=ax2,vmin=-20, vmax=20)#min_lon=klon1-offset,max_lon=klon1+offset,

    plt.savefig('/Users/severe2/Research/TORUS_Data/'+day+'/radar/Nexrad/plots/work'+str(i)+'_'+str(j)+'.png')

    return
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


##########################################
# Create new datasets for each platform ##
##########################################

NSSLMM_df, max_array, min_array = [], [], []
for MM in ['FFld','LIDR','Prb1','Prb2','WinS']:
    
    #uncomment to plot only a sub section of the time (I think) 
    #tstart, tend = dt.datetime(2019,5,21,0,50,0), dt.datetime(2019,5,21,1,0,0)
    #tstart, tend = None, None
    
    try: 
        MMfile='/Users/severe2/Research/TORUS_data/'+day+'/mesonets/NSSL/'+MM+'_'+day[2:]+'_QC_met.dat'
        print(MMfile) 
        if MM == 'FFld':
            FFld_df = read_nsslmm(MMfile,tstart,tend)
            FFld_df.name='FFld'
            NSSLMM_df.append(FFld_df)
            masked_df, unmasked_df = maskdata(p_var, FFld_df)
            max_val, min_val= masked_df.max(), masked_df.min()
        elif MM == 'LIDR':
            LIDR_df = read_nsslmm(MMfile,tstart,tend)
            LIDR_df.name='LIDR'
            NSSLMM_df.append(LIDR_df)
            masked_df, unmasked_df = maskdata(p_var, LIDR_df)
            max_val, min_val= masked_df.max(), masked_df.min()
        elif MM == 'Prb1':
            Prb1_df = read_nsslmm(MMfile,tstart,tend)
            Prb1_df.name='Prb1'
            NSSLMM_df.append(Prb1_df)
            masked_df, unmasked_df = maskdata(p_var, Prb1_df)
            max_val, min_val= masked_df.max(), masked_df.min()
        elif MM == 'Prb2':
            Prb2_df = read_nsslmm(MMfile,tstart,tend)
            Prb2_df.name='Prb2'
            NSSLMM_df.append(Prb2_df)
            masked_df, unmasked_df = maskdata(p_var, Prb2_df)
            max_val, min_val= masked_df.max(), masked_df.min()
        elif MM == 'WinS':
            WinS_df = read_nsslmm(MMfile,tstart,tend)
            WinS_df.name='WinS'
            NSSLMM_df.append(WinS_df)
            masked_df, unmasked_df = maskdata(p_var, WinS_df)
            max_val, min_val= masked_df.max(), masked_df.min()
        
        if MM == probe_of_interest:
            p_of_int=NSSLMM_df[-1]

    except: error_printing(e_test)
    
    try:
        max_array.append(max_val)
        min_array.append(min_val)
    except: error_printing(e_test)

    #max_array.append(max_val)
    #min_array.append(min_val)
#max an min values that will plot
#if p_var =='Thetae':
#    globalamin, globalamax=325, 360
#elif p_var =='Thetav':
#    globalamin, globalamax=295, 315

#Alternate global variable max/min for plotting purposes
globalamax = np.max(max_array)
globalamin = np.min(min_array)
print(p_of_int)
# ************************************

#print(Prb1_df)
#for i in len(Prb1_df):
#    Prb1_df['Radar_ID']=nearest_WSR(
test=det_nearest_WSR(Prb1_df)
r_ofintrest=test.Radar_ID.unique()

timeranges_each_r=pd.DataFrame()
i=0
for r in r_ofintrest:
    i=i+1
    print(r)
    t=test.loc[test.Radar_ID==r, ['datetime']].rename(columns={'datetime':r})
    ts,te=t.min(),t.max()
    timeranges_each_r=pd.concat([timeranges_each_r,t],axis=1)
    print("start ",ts[r])
    print("end ",te[r])
    print("***")
    radar_list=get_WSR_from_AWS(ts[r],te[r],r,temploc)
    j=0
    for radar in radar_list:
        j=j+1
        #print(radar.info())
        #testing_plots(radar,i,j,r_only,globalamin,globalamax,p_var,e_test)
        #plot_with_WSR(radar,i,j)

        ppiplot(r_only, radar_list,radar, filesys, day, globalamin,globalamax, p_var, p_of_int, e_test,i,j)
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
