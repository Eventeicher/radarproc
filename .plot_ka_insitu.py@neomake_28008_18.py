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
from read_platforms import read_unlmm
from read_platforms import read_nsslmm
import read_platforms
import sys, traceback

##########################################################################
#Variables
##########
day = '20190517' #'YYYYMMDD'
radar_plotting= True #would you like to plot radar images? (set to False to only plot timeseries)
r_only=False #set to true if you want the radar plots only 
p_var = "Thetav" #which var to plot (current options; Thetae, Thetav)
filesys='/Users/severe2/Research/'

#Troubleshooting y/n
####################
#Set to True to print statements when entering/leaving definitions (helpful to understanding whats going on in the code but prints alot of lines to your terminal)
print_long = False # True/False variable

#Set to False to not print out the errors that tripped 'try' statements: set to True only while troubleshooting 
  ##there will be 'valid' errors caused by looping through platforms etc.... hence needlessly confusing unless troubleshooting
e_test = False #True/False variable

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
#################################################
# Defintions that are not used (at the moment) ##
#################################################
def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)
    return
#* * * * * 
def resize_colorbar(event):
    plt.draw()
    posn = ax.get_position()
    cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0,0.04, posn.height])
    return
#* * * * * 
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
#* * * * * 
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return idx


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
    ''' Determine if a given radar is deployed and if so assign the correct values to it. 
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
        kadep, k_lat, k_lon, head_ing, rhi_b, rhi_e = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    
    return kadep, k_lat, k_lon, head_ing, rhi_b, rhi_e

# * * * * * * * 
def getRadarPositionData(c_time, r_dep):
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
def maskdata(p_var, NSSLMM):
    p_um = NSSLMM[p_var].values
    p_m = np.ma.masked_where(NSSLMM['qc_flag'].values>0, p_um)
    
    return p_m, p_um

#* * * * * 
def getLocation(lat1, lon1, brng, distancekm):
    lat1 = lat1 * np.pi / 180.0
    lon1 = lon1 * np.pi / 180.0
    
    R = 6378.1 #earth radius
    #R = ~ 3959 MilesR = 3959
    bearing = (brng / 90.)* np.pi / 2.

    lat2 = np.arcsin(np.sin(lat1) * np.cos(distancekm/R) + np.cos(lat1) * np.sin(distancekm/R) * np.cos(bearing))
    lon2 = lon1 + np.arctan2(np.sin(bearing)*np.sin(distancekm/R)*np.cos(lat1),np.cos(distancekm/R)-np.sin(lat1)*np.sin(lat2))
    lon2 = 180.0 * lon2 / np.pi
    lat2 = 180.0 * lat2 / np.pi
    
    return lat2, lon2

#* * * * * 
def createCircleAroundWithRadius(lat, lon, radiuskm,sectorstart,sectorfinish,heading):
    #    ring = ogr.Geometry(ogr.wkbLinearRing)
    latArray = []
    lonArray = []
    for brng in range(int(heading+sectorstart),int(heading+sectorfinish)): #degrees of sector
        lat2, lon2 = getLocation(lat,lon,brng,radiuskm)
        latArray.append(lat2)
        lonArray.append(lon2)
    
    return lonArray,latArray

# * * * * * * 
def rhi_spokes_rings(rhib,rhie,head,klat,klon,radar,display,ymin,ymax,xmin,xmax):
    if print_long== True: print('made it into rhi_spokes_rings')
    #produce spoke and ring
    for j in range(int(rhib),int(rhie)+1,10):
        ang= head + j
        if ang > 360.:
            ang=int(ang-360.)
        A,B = createCircleAroundWithRadius(klat,klon,(radar.range['data'][-1]-500.)/1000., rhib,rhie+1,head) #this plots a circle that connects the spokes
        C,D = getLocation(klat, klon, ang, (radar.range['data'][-1]-500.)/1000.)
        display.plot_line_geo(A,B,marker=None,color='grey',linewidth=.25) #this plots a circle that connects the spokes
        display.plot_line_geo([klon,D], [klat,C], marker=None, color='k', linewidth=0.5, linestyle=":")
        #if np.logical_and(C>ymin,np.logical_and(C<ymax,np.logical_and(D>xmin,D<xmax))):
            #d1=plt.text(D, C, str(ang),horizontalalignment='center',transform=ccrs.PlateCarree(),fontsize=10,zorder=9,path_effects=([PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')])
    if print_long== True: print('made it through rhi_spokes_rings')
    
    return

# * * * * * *
def legend_maker(p_name,m_style,m_color,leg_str,legend_elements):
    if print_long== True: print('made it into legend_maker')
    if (p_name in ['Ka1','Ka2']): #the legend entries for the KA radars
        legend_entry=Line2D([], [], marker=m_style, markeredgecolor='black',markeredgewidth=3,label=leg_str,markerfacecolor=m_color, markersize=26)
    else:
        legend_entry=Line2D([], [], marker=m_style, markeredgecolor=m_color,markeredgewidth=3,label=leg_str, markersize=26,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')])
    
    legend_elements=np.append(legend_elements,legend_entry)
    if print_long== True: print('made it through legend maker')
    
    return legend_elements

# * * * * * *
def platform_attr(p_file,legend_elements,radar_m=False,r_s=False):
    if print_long== True: print('Made it into platform_attr')
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
        p_name = p_file.name

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
        if print_long== True: print('The platform was not deployed at this time')
        p_deploy = False
        error_printing(e_test)

    return p_sub, lon_sub, lat_sub, U_sub, V_sub, p_deploy

# * * * * ** *    
def platform_plot(file,radartime,ax_n,color,m_color,p_var,e_test,labelbias=(0,0)):
    if print_long== True: print('made it into platform_plot')

    #grab the subset of data of +- interval around radar scan
    p_sub, lon_sub, lat_sub, U_sub, V_sub, p_deploy = grab_platform_subset(file, radartime, p_var)
    
    if p_deploy == True:
        #plot data for reflectivity plot
        plot_colourline(lon_sub,lat_sub,p_sub,cmocean.cm.curl,ax_n,ccrs.PlateCarree(),color=color,amin=globalamin,amax=globalamax)
  
        #d2=plt.text(lon_sub[int(len(lon_sub)/2)]+labelbias[0], lat_sub[int(len(lat_sub)/2)]+labelbias[1], file[58:62],#transform=ccrs.PlateCarree(), fontsize=20, zorder=9, path_effects=[PathEffects.withStroke(linewidth=4,foreground=color)])
        ax_n.plot(lon_sub[int(len(lon_sub)/2)], lat_sub[int(len(lat_sub)/2)], transform=ccrs.PlateCarree(), marker='1',markersize=20, markeredgewidth='3',color=m_color,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')], zorder=10)
  
        ax_n.plot(lon_sub[int(len(lon_sub))-1], lat_sub[int(len(lat_sub))-1], transform=ccrs.PlateCarree(), marker='.',markersize=10, markeredgewidth='3',color='k',zorder=9) 
  
        #plot windbarbs
        try:
            stationplot = metpy.plots.StationPlot(ax_n, lon_sub[::30], lat_sub[::30], clip_on=True, transform=ccrs.PlateCarree())
            stationplot.plot_barb(U_sub[::30], V_sub[::30],length=7)
        except: error_printing(e_test)

    if print_long== True: print('made it through platform_plot')
    
    return



##  #  #  #  #  #  #  #  #  #  #  #  #  # # # #  #  #  #  #  #  #  #  #  #  # #  #  #  # # #  #  #  #  #  #  #  #  # #  #
##  #  #  #  #  #  #  #  #  #  #  #  #  # # # #  #  #  #  #  #  #  #  #  #  # #  #  #  # # #  #  #  #  #  #  #  #  # #  #
def ppiplot(r_only,radar,swp_id,azimuth,currentscantime,klat,klon,klat1,klon1,klat2,klon2,rhib1,rhie1,rhib2,rhie2,head1,head2,globalamin,globalamax,p_var,e_test):
    if print_long== True: print('made it into ppiplot')
    
    SMALL_SIZE, MS_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 25, 30, 35, 50
    plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=MS_SIZE)       # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    ## Set up plot size
    if r_only == True:
        fig = plt.figure(figsize=(40,150), facecolor='white')
    else:
        fig = plt.figure(figsize=(30,20),facecolor='white')
    display = pyart.graph.RadarMapDisplay(radar)

    ## Plot instuments on subplots
    radar_subplots('refl',fig,display,klat,klon,klat1,klon1,klat2,klon2,rhib1,rhie1,rhib2,rhie2,head1,head2,radar,swp_id,currentscantime,p_var,e_test)
    radar_subplots('vel',fig,display,klat,klon,klat1,klon1,klat2,klon2,rhib1,rhie1,rhib2,rhie2,head1,head2,radar,swp_id,currentscantime,p_var,e_test)
    if r_only == False:
        time_series(filesys,day,fig,currentscantime,globalamin,globalamax,p_var,radar_sub=True)

    ## Plot colorbars
    if p_var == "Thetae":
        c_lab= "Equivalent Potential Temp [K]"
    elif p_var == "Thetav":
        c_lab= "Virtual Potential Temp [K]"
    cbar_ax = plt.axes([.5265,.5356,.014, 0.405])#left, bottom, width, height
    cbar=plt.colorbar(CS3,cax=cbar_ax, orientation='vertical', label=c_lab, ticks=MaxNLocator(integer=True))#,ticks=np.arange(globalamin,globalamax+1,2))
   
    ## Plot title
    title = thefile.split("/")[-1][10:13] + ' '+str(np.around(azimuth,2))+r'$^{\circ}$ PPI ' +currentscantime.strftime("%m/%d/%y %H:%M")+ ' UTC'
    plt.suptitle(title)
    
    ## Finish plot 
    plt.tight_layout()
    plt.savefig('/Users/severe2/Research/TORUS_Data/'+day+'/mesonets/plots/'+thefile.split("/")[-1][10:13]+'_'+currentscantime.strftime('%m%d%H%M')+'_'+p_var+'.png')
    plt.close()

    if print_long== True: print('made it through ppiplot')
    print('Done Plotting')
    print('******************************************************************************************************')
    
    return 

# * * * * * *  *
def radar_subplots(mom,fig,display,klat,klon,klat1,klon1,klat2,klon2,rhib1,rhie1,rhib2,rhie2,head1,head2,radar,swp_id,currentscantime,p_var,e_test):
    if print_long== True: print('made it into radar_subplots')
    
    ## SET UP PLOTTING CONTROLS
    NSSLmm, NEBmm, UASd = True, False, False #which other platforms do you wish to plot 
    rhi_ring= True #do you want the rhi spokes
    country_roads, hwys, county_lines, state_lines = False, True, False, False #background plot features
    legend_elements=[]
    ###########################

    ## SET UP VARS FOR EACH RADAR MOMENTS
    if mom == 'refl':
        pos, field=221, 'refl_fix'
        c_scale, c_label='pyart_HomeyerRainbow','Radar Reflectivity [dbz]'
        vminb, vmaxb= -30.,30.
        p_title, leg ='Reflectivity', True
    elif mom == 'vel':
        pos, field=222, 'vel_fix'
        c_scale, c_label='pyart_balance','Velocity [m/s]'
        vminb, vmaxb= -40.,40.
        p_title, leg ='Radial Velocity', False
    
    ## Bounding box, x km away from the radar in all directions 
    xmin, xmax = getLocation(klat,klon,270,21)[1], getLocation(klat,klon,90,21)[1]
    ymin, ymax = getLocation(klat,klon,180,21)[0], getLocation(klat,klon,0,21)[0]
    
    ## SET UP SUBPLOTS   
    ax_n=fig.add_subplot(pos,projection=display.grid_projection)
    ax_n.plot(klon,klat,transform=display.grid_projection)
    display.plot_ppi_map(field, swp_id,title_flag=False, cmap=c_scale, ax=ax_n, vmin=vminb, vmax=vmaxb, colorbar_label= c_label, min_lon=xmin, max_lon=xmax, min_lat=ymin, max_lat=ymax, embelish=False) 
    ax_n.text(.5,-.065,p_title,transform=ax_n.transAxes,horizontalalignment='center',fontsize=40) #the radar subplot titles
        
    ## PLOT RHI SPOKES
    if rhi_ring == True:
        try:
            rhi_spokes_rings(rhib1,rhie1,head1,klat1,klon1,radar,display,ymin,ymax,xmin,xmax) #plot the spokes for radar1
        except: error_printing(e_test)
        try:
            rhi_spokes_rings(rhib2,rhie2,head2,klat2,klon2,radar,display,ymin,ymax,xmin,xmax) #plot the spokes for radar2
        except: error_printing(e_test)

    ## PLOT RADAR MARKERS
    try:
        if np.logical_and(klat2>ymin,np.logical_and(klat2<ymax,np.logical_and(klon2>xmin,klon2<xmax))):
            m_style,m_color,l_color,leg_str,legend_elements=platform_attr('Ka2',legend_elements,radar_m=True, r_s=True)
            ax_n.plot(klon2,klat2,marker=m_style, transform=ccrs.PlateCarree(),color=m_color,markersize=18,markeredgewidth=5,path_effects=[PathEffects.withStroke(linewidth=15,foreground='k')],
                    zorder=10)
            #d2=ax_n.text(klon2+0.006, klat2-0.011,'Ka2',transform=ccrs.PlateCarree(), zorder=10, path_effects=[PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')]) #textlabel on plot 
    except: error_printing(e_test)
    try:
        if np.logical_and(klat1>ymin,np.logical_and(klat1<ymax,np.logical_and(klon1>xmin,klon1<xmax))):
            m_style,m_color,l_color,leg_str,legend_elements=platform_attr('Ka1',legend_elements,radar_m=True, r_s=True)
            ax_n.plot(klon1,klat1,marker=m_style, transform=ccrs.PlateCarree(),color=m_color,markersize=18,markeredgewidth=5,path_effects=[PathEffects.withStroke(linewidth=15,foreground='k')],
                    zorder=10)
            #d1=ax_n.text(klon1+0.006, klat1+0.005, 'Ka1',transform=ccrs.PlateCarree(), zorder=10, path_effects=[PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')])
    except: error_printing(e_test)
    
    ## PLOT OTHER PLATFORMS             
    if NSSLmm == True:
        for NSSLMM in NSSLMM_df:
            if print_long== True: print(NSSLMM.name)
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
    
    ## PLOT BACKGROUND FEATURES
    if country_roads == True:
        ox.config(log_file=True, log_console=True, use_cache=True)
        G = ox.graph_from_bbox(ymax,ymin,xmax,xmin)
        ox.save_load.save_graph_shapefile(G, filename='tmp'+str(0), folder=filesys+'Radar_Processing/TORUS19/roads/', encoding='utf-8')
        fname = filesys+'Radar_Processing/TORUS19/roads/tmp'+str(0)+'/edges/edges.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree(), edgecolor='gray', linewidth=0.5)
        ax_n.add_feature(shape_feature, facecolor='none')
        shutil.rmtree(filesys+'Radar_Processing/TORUS19/roads/tmp'+str(0)+'/')
    if hwys == True:
        fname = filesys+'Radar_Processing/TORUS19/roads/GPhighways.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree(), edgecolor='grey')#edgecolor='black')
        ax_n.add_feature(shape_feature, facecolor='none')
    if county_lines == True:
        fname = filesys+'Radar_Processing/TORUS19/roads/cb_2017_us_county_5m.shp'
        shape_feature = ShapelyFeature(Reader(fname).geometries(),ccrs.PlateCarree(), edgecolor='gray')
        ax_n.add_feature(shape_feature, facecolor='none', linewidth=1.5, linestyle="--")
    if state_lines == True:
        states_provinces = cartopy.feature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
        ax_n.add_feature(states_provinces, edgecolor='black', linewidth=2)
     
    if print_long== True: print('Made it through radar_subplots')
    
    return

# * * * * * * * 
def time_series(filesys,day,fig,currentscantime,globalamin,globalamax,p_var,radar_sub=False):
    if print_long== True: print('Made it into time_series')
    if radar_sub == True:
        ax_n, var1='ax3', 212
        vline=currentscantime 
        ax_n=fig.add_subplot(var1)
    
    else: #will only plot the time series (not the radar subplots)
        ax_n='ax'
        fig, ax_n= plt.subplots()
    
    #NSSLMM_files= sorted(glob.glob(filesys+'TORUS_Data/'+day+'/mesonets/NSSL/*.nc'))
    for NSSLMM in NSSLMM_df:
        if print_long== True: print(str(NSSLMM.name))
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
    
    if print_long== True: print('Made it through time_series')
    
    return


###############################################################################################
##########################################
# Create new datasets for each platform ##
##########################################
#crop the start/end time to a time you specify (or set to None to get the full time range)
tstart = dt.datetime(int(day[0:4]),int(day[4:6]),int(day[6:8]),15,0,0)
tend= None

NSSLMM_df, UNLMM_df, max_array, min_array = [], [], [], []
for MM in ['FFld','LIDR','Prb1','Prb2','WinS','CoMeT1','CoMeT2','CoMeT3']:
   
    #load the NSSL mesonets 
    if (MM in ['FFld','LIDR','Prb1','Prb2','WinS']):
        try: 
            MMfile='/Users/severe2/Research/TORUS_data/'+day+'/mesonets/NSSL/'+MM+'_'+day[2:]+'_QC_met.dat'
        
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
                #max_val, min_val= WinS_df[p_var].max(), WinS_df[p_var].min()
        except: error_printing(e_test)

    #load the UNL MM files 
    elif (MM in ['CoMeT1', 'CoMeT2', 'CoMeT3']):
        try: 
            MMfile='/Volumes/Samsung_T5/Research/TORUS_Data/'+day+'/mesonets/UNL/UNL.'+MM+'.'+day+'*'
            
            if MM == 'CoMeT1':
                CoMeT1_df, co1 = read_unlmm(MMfile,tstart,tend)
                CoMeT1_df.name='CoMeT1'
                UNLMM_df.append(CoMeT1_df)
                #masked_df, unmasked_df = maskdata(p_var, CoMeT1_df)
                #max_val, min_val= unmasked_df.max(), unmasked_df.min()
            elif MM == 'CoMeT2':
                CoMeT2_df,co2 = read_unlmm(MMfile,tstart,tend)
                CoMeT2_df.name='CoMeT2'
                UNLMM_df.append(CoMeT2_df)
                #masked_df, unmasked_df = maskdata(p_var, CoMeT2_df)
                #max_val, min_val= unmasked_df.max(), unmasked_df.min()
            elif MM == 'CoMeT3':
                CoMeT3_df , co3= read_unlmm(MMfile,tstart,tend)
                CoMeT3_df.name='CoMeT3'
                UNLMM_df.append(CoMeT3_df)
                #masked_df, unmasked_df = maskdata(p_var, CoMeT3_df)
                #max_val, min_val= unmasked_df.max(), unmasked_df.min()
        except: error_printing(e_test)

    max_array.append(max_val)
    min_array.append(min_val)
print(co1)

#max an min values that will plot
#if p_var =='Thetae':
#    globalamin, globalamax=325, 360
#elif p_var =='Thetav':
#    globalamin, globalamax=295, 315

#Alternate global variable max/min for plotting purposes
globalamax = np.max(max_array)
globalamin = np.min(min_array)

print(UNLMM_df)
#######################################################################################################
#################
# Create Plots ##
#################
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
                        ka1dep, klat1, klon1, head1, rhib1, rhie1= det_radar_deps(currentscantime,'ka1',e_test)
                    except: error_printing(e_test)
                    try:
                        ka2dep, klat2, klon2, head2, rhib2, rhie2= det_radar_deps(currentscantime,'ka2',e_test)
                    except: error_printing(e_test)

                    #determine which radar is the main plotting radar
                    if thefile.split("/")[-1][10:13] == 'Ka1':
                        klat, klon, head = klat1, klon1, head1
                    elif thefile.split("/")[-1][10:13] == 'Ka2':
                        klat, klon, head = klat2, klon2, head2
                    else:
                        print('What radar are we using?')

                    #proceed to plot the radar
                    ppiplot(r_only,radar,swp_id,azimuth,currentscantime,klat,klon,klat1,klon1,klat2,klon2,rhib1,rhie1,rhib2,rhie2,head1,head2,globalamin,globalamax,p_var,e_test)

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
#Only plotting timeseries
    print("no radar plotted")
    time_series(filesys,day,fig,p_var,radar_sub=False)
    #print(str(filesys+'TORUS_Data/'+day+'/mesonets/NSSL/*.nc'))

    fig.savefig('test2.png')
    plt.close()
