#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import needed modules
######################
import matplotlib
#  matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.patheffects as PathEffects
from matplotlib import ticker
from matplotlib.ticker import (LinearLocator, FixedLocator, MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.cm as cmx
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler
import datetime as dt
from datetime import datetime
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import ShapelyFeature, NaturalEarthFeature
from cartopy.io.shapereader import Reader
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
import metpy
from metpy.plots import StationPlot, ctables
from metpy.units import units
import metpy.calc as mpcalc
#  import modin.pandas as pd
import pandas as pd
from pandas.plotting import register_matplotlib_converters
import numpy as np
import xarray as xr
#from netCDF4 import num2date
import logging
import os
import os.path
from joblib import Memory, Parallel, delayed
from scipy import ndimage, interpolate
from collections import namedtuple
import cmocean
import gc
import glob
import nexradaws
import pyart
import shutil
import sys
import traceback
import cftime # for num2pydate
import pprint
import math

pp = pprint.PrettyPrinter(indent=4)

from radar_defns import radar_fields_prep, fix_NOXP_order 
register_matplotlib_converters()
logging.basicConfig(filename='this_is.log')
log = logging.getLogger(__name__)

################################################################################################
##################
# TROUBLE SHOOTING
##################
def error_printing(e_test):
    ''' Basically I got sick of removing/replacing the try statements while troubleshooting '''
    if e_test == True:
        e = sys.exc_info()
        traceback.print_tb(e[-1])
        #print( "<p>Error: %s</p>" % e )
        print("Error:", e[:-2], '\n')

def timer(start, end, total_runtime=False):
   hours, rem = divmod(end-start, 3600)
   minutes, seconds = divmod(rem, 60)
   if total_runtime == False:
       print("Plot took {:0>2} hrs, {:0>2} min, {:05.2f} sec to complete".format(int(hours),int(minutes), seconds))
   if total_runtime == True:
       print("\nIt took {:0>2} hrs, {:0>2} min, {:05.2f} to complete all of the plots\n".format(int(hours),int(minutes), seconds))

###########
# Data Prep
###########
def time_range(config):
    #if you specify a start or end time it will be assigned here otherwise will be set to none (full data set)
    try:  Tstart = config.tstart
    except: Tstart = None

    try: Tend = config.tend
    except: Tend = None

    return Tstart, Tend

def time_in_range(start, end, x):
    """Return true if x is in the range [start, end]"""
    #  if end == None: end=datetime.utcnow()
    if end == None: end=datetime.max
    if start == None: start=datetime.min
    if start <= end: return start <= x <= end
    else: return start <= x or x <= end

def pform_names(Type):
    ''' Returns lists of platform names for various platform types; Saves typing in the end
        ---
        INPUT: Type [str]: valid options include 'ALL','RADAR', 'TInsitu','UNL','NSSL','KA'
        Output: P_list: list of pnames containing the names of the requested type
    '''
    if Type == 'ALL': P_list = ['WTx_M','OK_M','IA_M','KS_M','ASOS','FFld','LIDR','Prb1','Prb2','WinS','CoMeT1','CoMeT2','CoMeT3','UAS', 'Ka1','Ka2','NOXP','WSR88D']
    elif Type == "STN_I": P_list= ['WTx_M','OK_M','IA_M','KS_M', 'ASOS']
    elif Type == 'MESO': P_list = ['WTx_M','OK_M','IA_M', 'KS_M']
    elif Type == 'TInsitu': P_list = ['FFld','LIDR','Prb1','Prb2','WinS','CoMeT1','CoMeT2','CoMeT3']#'UAS']
    elif Type == 'UNL': P_list = ['CoMeT1','CoMeT2','CoMeT3']
    elif Type == 'NSSL': P_list = ['FFld','LIDR','Prb1','Prb2','WinS']
    elif Type == 'RADAR': P_list = ['Ka1','Ka2','NOXP','WSR88D']
    elif Type == 'KA': P_list = ['Ka1','Ka2']
    else: print('Please enter a valid list name')
    return P_list

# *******************
def pform_attr(pname):
    ''' INPUTS: pname: name of the platform
        ----
        RETURNS: legend_entry: The specs regarding the legend presentation for each pform (will be appended latter to make the legend)
                 m_color,m_style etc : the specifications regarding platform presentation on the plots (markershape, color, etc)
    '''
    ##assign the atributes for each platform
    if pname in pform_names('TInsitu'): 
        #  marker_style, marker_size = '$\U0001F766$', 20 #u20B9 u03F0 u03A8 u273D 204E 0359 002A 20BA 2723
        marker_style, marker_size = '$\u2724$', 10 #u20B9 u03F0 u03A8 u273D 204E 0359 002A 20BA 2723
        #  marker_style, marker_size = '$\u2724$', 20 #u20B9 u03F0 u03A8 u273D 204E 0359 002A 20BA 2723
        if pname == "Prb1": marker_color, line_color, legend_str= 'xkcd:lightblue', 'steelblue', 'Prb1'
        elif pname == "Prb2": marker_color, line_color, legend_str= 'xkcd:watermelon', 'xkcd:dusty red', 'Prb2'
        elif pname == "FFld": marker_color, line_color, legend_str= 'xkcd:bubblegum pink', 'xkcd:pig pink', 'FFld'
        elif pname == "LIDR": marker_color, line_color, legend_str= 'xkcd:pastel purple', 'xkcd:light plum', 'LIDR'
        elif pname == "WinS": marker_color, line_color, legend_str= 'xkcd:peach', 'xkcd:dark peach', 'WinS'
        elif pname == "CoMeT1": marker_color, line_color, legend_str= 'springgreen', 'mediumseagreen','CoMeT1'
        elif pname == "CoMeT2": marker_color, line_color, legend_str= 'yellow', 'gold', 'CoMeT2'
        elif pname == "CoMeT3": marker_color, line_color, legend_str= 'peachpuff', 'sandybrown', 'CoMeT3'

    if pname in pform_names('STN_I'):
        marker_color, line_color = 'black', 'black'
        if pname in pform_names('MESO'):
            marker_size = 13 #U+1278, #u20a9
            #  marker_size = 23 #U+1278, #u20a9
            if pname == "WTx_M": marker_style,  legend_str= '$\u0166$', 'WTxM'
            elif pname == "OK_M": marker_style, legend_str= '$\u00C5$', 'OKM' ##00F0 2644 26B7  26A8 0516 2114
            elif pname == "IA_M": marker_style, legend_str= '$\uABA0$', 'IAM'
            elif pname == "KS_M": marker_style, legend_str= '$\u0199$', 'KSM'#0199 A741 0198 01E9 03CF ##0198 0199
        else:
            marker_size = 15
            #  marker_size = 25
            if pname == "ASOS": marker_style, legend_str= '$\u0466$', 'ASOS' #0466 212B 00C4 20B3 00C6 0102 01E0 01CD 01DE

    if pname in pform_names('RADAR'):
        if pname in pform_names('KA'): ##2643
            marker_style, marker_size = '$\u264E$', 20
            #  marker_style, marker_size = '$\u264E$', 35
            #  marker_style, marker_size = '$\u14C7$', 35
            if pname == 'Ka2': marker_color, line_color, legend_str= 'xkcd:crimson', 'xkcd:crimson', 'Ka2'
            elif pname == 'Ka1': marker_color, line_color, legend_str= 'mediumseagreen', 'mediumseagreen', 'Ka1'
        else:
            line_color, marker_size = 'black', 20
            #  line_color, marker_size = 'black', 35
            if pname == 'NOXP': marker_style, marker_color, legend_str= '$\u2116$', 'yellow', "NOXP"
            #  if pname == 'NOXP': marker_style, marker_color, legend_str= '$\u306E$', 'green', "NOXP"
            elif pname == 'WSR88D':  
                marker_style, legend_str= '$\U0001D6C0$', "WSR88D"
                #  marker_color = (cycler(color=['r', 'g', 'b', 'y']))
                marker_color = 'white' 
                  
       ##create the legend elements that will be appended to the legend array
    if pname in pform_names('KA'): #the legend entries for the KA radars
        legend_entry = Line2D([], [], marker=marker_style, markeredgecolor='black', markeredgewidth=3, label=legend_str, markerfacecolor=marker_color, markersize=26)
    elif pname in ['WTx_M','OK_M','IA_M','KS_M', 'ASOS']:
        legend_entry = Line2D([], [], marker=marker_style, markeredgecolor=marker_color, label=legend_str, markersize=26)
    else:
        legend_entry = Line2D([], [], marker=marker_style, markeredgecolor=marker_color, markeredgewidth=3, label=legend_str, markersize=26, 
                              path_effects=[PathEffects.withStroke(linewidth=12, foreground='k')])
    return marker_style, marker_color, marker_size, line_color, legend_str, legend_entry

################################################################################################
#  * * * * *
def Add_to_DATA(config, day, DType, Data, subset_pnames, MR_file=None, swp=None):
    if config.print_long == True: print('Made it into Add_to_DATA')
    if DType == 'TInsitu':
        #for these pforms mainly interested if we have any valid data for the day not worried about spatial/temporal subset yet
        #  To save computing time this data will only be read in once and anyfurther subsetting will be done later
        for pname in pform_names('TInsitu'):
            # if you do want to read in data
            if config.Read_Data['Meso'] == True:
                #for each platform test to see if we have data for that day
                data_avail = Platform.test_data(config, day, pname, scan_time= None)
                if data_avail == True:
                    subset_pnames.append(pname) #append the pname to the subset_pnames list
                    #  load data for the pform (aka initialize an object of the appropriate class); place in dict with key of pname
                    pform_data= Torus_Insitu(config, day, pname)
                    Data.update({pname: pform_data})
                    if config.print_long == True: print("Data can be read in for platform %s" %(pname))
                else:
                    if config.print_long == True: print("No data available to be read in for platform %s" %(pname))
            elif config.Read_Data['Meso'] == False:
                    if config.print_long == True: print("Did not attempt to read in data for platform %s" %(pname))

    # * * *
    elif DType == 'STN_I':
        #read in the site location for mesonet datasets (W Tx, OK, IA) .
        #  Will only be called once since the position of the sites wont change for each image
        for pname in pform_names('STN_I'):
            # if you do want to read in data
            if config.Read_Data['Meso'] == True:
                data_avail = Platform.test_data(config, day, pname, scan_time=None)
                if data_avail == True:
                    #  dont repeatedly append to the list if the platform is already included
                    if pname in subset_pnames: pass
                    else: subset_pnames.append(pname) #append the pname to the subset_pnames list
                    #  load loc data for the sites in plotting domain
                    pform_data=Stationary_Insitu(config, day, pname)
                    Data.update({pname: pform_data})
                    if config.print_long == True: print("Data can be read in for platform %s" %(pname))
                else:
                    if config.print_long == True: print("No data available to be read in for platform %s" %(pname))
            elif config.Read_Data['Meso'] == False:
                    if config.print_long == True: print("Did not attempt to read in data for platform %s" %(pname))

    # * * *
    elif DType == 'RADAR':
        #this can change for each image (aka from diff scans of the plotting radar) so will be updated for every plot
        #  ie. are the radars deployed at a given time which (if any) WSR88D fall within plotting domain etc
        # * * *
        ## Initilize the Main plotting radar object (which will establish scantime and contain the radar data file etc)
        # Det Main plotting radar (aka which radar is associated with MR_file)
        if MR_file != None:
            if MR_file.metadata['instrument_name'] == 'TTUKa-1': MR_name = 'Ka1'
            elif MR_file.metadata['instrument_name'] == 'TTUKa-2': MR_name = 'Ka2'
            elif MR_file.metadata['instrument_name'] == 'NOXPRVP': MR_name = 'NOXP'
            elif MR_file.metadata['original_container'] == 'NEXRAD Level II': MR_name = 'WSR88D'
            else: print('What radar are you trying to plot? MR_name = %' %(MR_file.metadata['instrument_name']))

            print("Reading in radar data for plotting from %s" %(MR_name))
            pform_data= Radar(config, day, MR_name, Rfile=MR_file, Swp=swp, Plotting_Radar=True)
            Data.update({'P_Radar': pform_data})

        # * * *
        ## Initilize the other radar objects that do not contain radar data to be plotted
        #  (aka could include marker locations for the other platforms etc)
        for pname in pform_names('RADAR'):
            # if you do want to read in data
            if config.Read_Data['Radar'] == True:
                data_avail = Platform.test_data(config, day, pname, scan_time=Data['P_Radar'].Scan_time)
                if data_avail == True:
                    #  dont repeatedly append to the list if the platform is already included
                    if pname in subset_pnames: pass
                    else: subset_pnames.append(pname) #append the pname to the subset_pnames list
                    #  load loc data (and in this case rhi angles if applicable)
                    pform_data= Radar(config, day, pname, Data)
                    Data.update({pname:pform_data})
                    if config.print_long == True: print("Data can be read in for platform %s" %(pname))
                else:
                    if config.print_long == True: print("No data available to be read in for platform %s" %(pname))
            elif config.Read_Data['Radar'] == False:
                    if config.print_long == True: print("Did not attempt to read in data for platform %s" %(pname))
    # * * *
    ##Uncomment to check yourself
    # ***************************
    #  print(subset_pnames)
    #  print(Data)
    #  print(dir(Data['FFld']))
    #  print(vars(Data['FFld']))
    #  print(Data['FFld'].m_color)
    #  print(Data[p_var].global_max)

    ##A few useful built in functions
    # *******************************
    #  dir(object): return all the properties and methods (including methods that are built in by default) of an object
    #  vars(object): returns all the values of the attributes for a given object
    #  getattr(abject,attribute): get the attibute value of the object, if object doesnt have the attr can return a default value
    #  hasattr(object, attribute): like getattr but returns a true/false bool
    #  isinstance, issubclass
    #  there are more: check google for built in class functionalitly
    if config.print_long == True: print('Made it through Add_to_DATA')
    return Data, subset_pnames


################################################################################################
#**************
def read_TInsitu(config, day, pname, d_testing=False):
    ''' Reads data files provided on TORUS19 EOL site (Important that datafiles and filenames follow TORUS19 readme)
        ---
        INPUT: filename string (following readme conventions for each platform)
               Tstart,Tend: datetime objects
        OUTPUT: dataframe containing all data collected during desired times
    '''
    def get_timeave(df, var, start, end):
        df_sub = df.loc[(df['datetime'] >= start) & (df['datetime'] <= end)]
        pd.set_option("display.max_rows", None)

        #  print(df_sub[['datetime',var]])
        mean = df_sub[var].mean()
        return mean
    def get_timearround(df, var, centered_time, time_offset):
        lower_tbound = (centered_time-dt.timedelta(seconds=time_offset)) 
        upper_tbound = (centered_time+dt.timedelta(seconds=time_offset)) 
        df_sub = df.loc[(df['datetime'] >= lower_tbound) & (df['datetime'] <= upper_tbound)]
        return df_sub 
    def calc_slope(x):
        slope = np.polyfit(range(len(x)), x, 1)[0]
        return slope
    ####
    if pname in pform_names('UNL'):
        mmfile = glob.glob(config.g_TORUS_directory+day+'/data/mesonets/UNL/UNL.'+pname+'.*')

        # Test Data availability
        if d_testing == True:
            try:
                #if there is no files this will will cause the try statement to fail
                mmfile[0]
                return True
            except: return False
        #  if the defn was only called to it will exit at this point (saving computing time)
        # + + + + + + + + + + + + ++ + +

        print(mmfile)
        data_hold = [] #empty list to append to
        for i in range(len(mmfile)):
            ds = xr.open_dataset(mmfile[i])

            #convert form epoch time to utc datetime object
            timearray = np.array([dt.datetime.utcfromtimestamp(int(t)/1e9) for t in ds.time.values])
            U,V = mpcalc.wind_components(ds.wind_speed, ds.wind_dir)
            lats = np.array([40]) * units.degrees
            
            #create new xarray dataset to make plotting easier 
            dims = ['datetime']
            coords = { 'datetime': timearray }
            data_vars = {
                'lat': (dims, ds.lat.values, {'units':str(lats.units)}),
                'lon': (dims, ds.lon.values, {'units':str(lats.units)}),
                'Z_ASL': (dims, ds.alt.values, {'units':str(ds.alt.units)}),
                'Z_AGL': (dims, np.zeros_like(ds.alt.values), {'units':str(ds.alt.units)}),
                'tfast': (dims, ds.fast_temp.values, {'units':str(ds.fast_temp.units)}),
                'Dewpoint': (dims, ds.dewpoint.values, {'units':str(ds.dewpoint.units)}),
                #  'RH': (dims, ds.calc_corr_RH.values, {'units':str(ds.calc_corr_RH.units)}),
                'Pressure': (dims, ds.pressure.values, {'units':str(ds.pressure.units)}),
                'spd': (dims, ds.wind_speed.values, {'units':str(ds.wind_speed.units)}),
                'dir': (dims, ds.wind_dir.values, {'units':str(ds.wind_dir.units)}),
                'U': (dims, U.data, {'units':str(ds.wind_speed.units)}),
                'V': (dims, V.data, {'units':str(ds.wind_speed.units)}),
                'Theta': (dims, ds.theta.values, {'units':str(ds.theta.units)}),
                'Thetae': (dims, ds.theta_e.values, {'units':str(ds.theta_e.units)}),
                'Thetav': (dims, ds.theta_v.values, {'units':str(ds.theta_v.units)})
            }

            subds = xr.Dataset(data_vars, coords)

            #convert to pandas
            pd_unl = subds.to_dataframe()
            pd_unl.reset_index(inplace=True)

            # Subset the dataset to the desired 22:43
            time_start, time_end= time_range(config)
            if time_start is None: hstart = pd_unl['datetime'].iloc[0]
            else: hstart = time_start
            if time_end is None: hend = pd_unl['datetime'].iloc[-1]
            else: hend = time_end
            # Save only desired iop data
            data_u = pd_unl.loc[(pd_unl['datetime'] >= hstart) & (pd_unl['datetime'] <= hend)]
            data_hold.append(data_u)

        #convert the list holding the dataframes to one large dataframe
        data_unl = pd.concat(data_hold)

        if config.overlays['Colorline']['Pert'] == True:
            for var in config.Line_Plot: 
                base_times=config.Pert_times[day]
                var_pert_data=[]
                if base_times != None:
                    for p in base_times.keys(): 
                        if p == 'All':
                            tstrt= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times['All']['starth'], base_times['All']['startm'], 0)
                            tend= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times['All']['endh'], base_times['All']['endm'], 0)
                            base_state = get_timeave(data_unl, var, tstrt, tend) 
                        elif p == 'All except p1':
                            if pname != 'Prb1':
                                tstrt= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times['All except p1']['starth'], base_times['All except p1']['startm'], 0)
                                tend= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times['All except p1']['endh'], base_times['All except p1']['endm'], 0)
                                base_state = get_timeave(data_unl, var, tstrt, tend) 
                            else: 
                                base_state=None
                        elif pname == p:
                            tstrt= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times[p]['starth'], base_times[p]['startm'], 0)
                            tend= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times[p]['endh'], base_times[p]['endm'], 0)
                            base_state = get_timeave(data_unl, var, tstrt, tend) 
                        else: 
                            base_state=None

                        if base_state != None:
                            #  pd.set_option("display.max_columns", None)
                            #  print(data_nssl)
                            print(('Base: {}, Pform: {}').format(base_state, pname)) 
                            for i,j in data_unl.iterrows():
                                #  print(('Base: {}, orig: {}, diff: {}, time: {}').format(base_state, j[var], base_state-j[var], j['datetime']))
                                pert_data= j[var] - base_state 
                                var_pert_data.append(pert_data)
                        else: 
                            var_pert_data=np.full_like(data_unl[var], np.nan)
                else:
                    var_pert_data=np.full_like(data_unl[var], np.nan)

                data_unl.loc[:, var+str('_pert')]=var_pert_data    
        if config.lineplt_control['Deriv'] == True:
            for var in config.Line_Plot: 
                if config.overlays['Colorline']['Pert'] == True:
                    var= var+str('_pert')
                der_data=[]
                for i,j in data_unl.iterrows():
                    p_sub=get_timearround(data_unl, var, j['datetime'], 30)
                    slope=calc_slope(p_sub[var])
                    der_data.append(slope)
                data_unl.loc[:, var+str('_der')]=der_data    


        #drop all the columns that we will not use past this point (to save memory/computing time)
        data_unl = data_unl.drop(columns=['Z_ASL', 'Z_AGL', 'Theta', 'Dewpoint', 'Pressure'])#, 'RH'])
        #  print(data_unl.memory_usage())
        #  data_unl.set_index('datetime', inplace=True, drop=False)
        return data_unl, 'UNL'

    # * * *
    elif pname in pform_names('NSSL'):
        mmfile = glob.glob(config.g_TORUS_directory+day+'/data/mesonets/NSSL/'+pname+'_'+day[2:]+'_QC_met.dat')
        
        # Test Data availability
        if d_testing == True:
            try:
                #if there is no files this will will cause the try statement to fail
                mmfile[0]
                return True
            except: return False
        # + + + + + + + + + + + + ++ + +
        mmfile=mmfile[0]
        # Read NSSL file using column names from readme
        column_names = ['id','time','lat','lon','alt','tfast','tslow','rh','p','dir','spd','qc1','qc2','qc3','qc4']
        
        data = pd.read_csv(mmfile, header=0, delim_whitespace=True, names=column_names)
        data = data.drop_duplicates()

        # Find timedelta of hours since start of iop (IOP date taken from filename!)
        print(mmfile)
        tiop = dt.datetime(2019, int(mmfile[-15:-13]), int(mmfile[-13:-11]), 0, 0, 0)

        time_start, time_end= time_range(config)
        if time_start is None:
            hstart = tiop
            hstart_dec = hstart.hour + (hstart.minute/60) + (hstart.second/3600) #convert to decimal hours HH.HHH
        else:
            hstart = float((time_start - tiop).seconds)/3600
            hstart_dec = hstart

        if time_end is None:
            hend = data['time'].iloc[-1]
        else:
            hend = (time_end - tiop)
            if hend >= dt.timedelta(days=1): hend = float((time_end-tiop).seconds)/3600 + 24.
            else:  hend = float((time_end-tiop).seconds)/3600

        # Save only desired iop data
        data_nssl, data_nssl_sub = pd.DataFrame(), pd.DataFrame()
        data_nssl_sub= data.loc[(data.loc[:,'time'] >= hstart_dec)]
        data_nssl= data_nssl_sub.loc[(data_nssl_sub.loc[:,'time'] <= hend)]
        # Convert time into timedeltas
        date = dt.datetime.strptime('2019-'+mmfile[-15:-13]+'-'+mmfile[-13:-11],'%Y-%m-%d')
        # Convert deltas into actual times
        data_nssl.loc[:,'datetime'] = pd.to_timedelta(data_nssl.loc[:,'time'], unit="h") + date

        ## Caclulate desired variables
        p, t = data_nssl['p'].values * units.hectopascal, data_nssl['tfast'].values * units.degC
        r_h = data_nssl['rh'].values/100
        data_nssl.loc[:, 'Theta'] = (mpcalc.potential_temperature(p, t)).m

        mixing = mpcalc.mixing_ratio_from_relative_humidity(p, t, r_h)
        data_nssl.loc[:, 'Thetav'] = (mpcalc.virtual_potential_temperature(p, t, mixing)).m

        td = mpcalc.dewpoint_from_relative_humidity(t, r_h)
        data_nssl.loc[:, 'Thetae'] = (mpcalc.equivalent_potential_temperature(p, t, td)).m

        Spd, dire = data_nssl['spd'].values * units('m/s') , data_nssl['dir'].values * units('degrees')
        u, v = mpcalc.wind_components(Spd, dire)
        data_nssl.loc[:, 'U'], data_nssl.loc[:, 'V'] = u.to('knot'), v.to('knot')

        if config.overlays['Colorline']['Pert'] == True:
            for var in config.Line_Plot:
                base_times=config.Pert_times[day]
                var_pert_data=[]
                #  print(vars(base_times))
                if base_times != None:
                    for p in base_times.keys(): 
                        if p == 'All':
                            tstrt= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times['All']['starth'], base_times['All']['startm'], 0)
                            tend= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times['All']['endh'], base_times['All']['endm'], 0)
                            base_state = get_timeave(data_nssl, var, tstrt, tend) 
                        elif p == 'All except p1':
                            if pname != 'Prb1':
                                tstrt= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times['All except p1']['starth'], base_times['All except p1']['startm'], 0)
                                tend= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times['All except p1']['endh'], base_times['All except p1']['endm'], 0)
                                base_state = get_timeave(data_nssl, var, tstrt, tend) 
                            else: 
                                base_state=None

                        elif pname == p:
                            tstrt= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times[p]['starth'], base_times[p]['startm'], 0)
                            tend= dt.datetime(int(day[0:4]), int(day[4:6]), int(day[6:8]), base_times[p]['endh'], base_times[p]['endm'], 0)
                            base_state = get_timeave(data_nssl, var, tstrt, tend) 
                        else: 
                            base_state=None

                        #  print(pname)
                        #  print(p)
                        #  print(base_state)
                        #  print('yyyyyy')
                        if base_state != None:
                            #  pd.set_option("display.max_columns", None)
                            #  print(data_nssl)
                            print(('Base: {}, Pform: {}').format(base_state, pname)) 
                            for i,j in data_nssl.iterrows():
                                #  print(('Base: {}, orig: {}, diff: {}, time: {}').format(base_state, j[var], j[var]-base_state, j['datetime']))
                                pert_data= j[var] - base_state 
                                var_pert_data.append(pert_data)
                        else: 
                            var_pert_data=np.full_like(data_nssl[var], np.nan)
                else: 
                    var_pert_data=np.full_like(data_nssl[var], np.nan)

                data_nssl.loc[:, var+str('_pert')]=var_pert_data    
        if config.lineplt_control['Deriv'] == True:
            for var in config.Line_Plot: 
                if config.overlays['Colorline']['Pert'] == True:
                    var= var+str('_pert')
                der_data=[]
                for i,j in data_nssl.iterrows():
                    p_sub=get_timearround(data_nssl, var, j['datetime'], 30)
                    slope=calc_slope(p_sub[var])
                    der_data.append(slope)
                data_nssl.loc[:, var+str('_der')]=der_data    


        q_list = config.overlays['NSSL']['qcflags'] 
        data_nssl.loc[:, 'all_qc_flags'] = data_nssl[q_list].sum(axis=1)
        #  data_nssl.set_index('datetime', inplace=True, drop=False)

        #drop all the columns that we will not use past this point (to save memory/computing time)
        data_nssl = data_nssl.drop(columns=['rh', 'p', 'time', 'alt', 'Theta', 'tslow'])
        #  print(data_nssl.memory_usage())
        return data_nssl, 'NSSL'

    # * * *
    elif pname == 'UAS':
        if config.print_long == True: print("no code for reading UAS yet")
        if d_testing == True: return False
        return 'UAS'

#**************
def read_Stationary(config, day, pname, d_testing=False):
    ''' Determine if a there are any sites from the stationary arrays that fall within the plot domain
        if so record locations and names of the sites
    '''
    if pname == 'WTx_M': file_name= 'West_TX_mesonets.csv'
    elif pname == 'OK_M': file_name= 'OKmeso.csv'
    elif pname == 'KS_M': file_name= 'KS_Meso.csv'
    elif pname == 'IA_M': file_name= 'IA_meso.csv'
    elif pname == 'ASOS': file_name= 'ASOS_stations.csv'

    stnry_df = pd.read_csv(config.g_csv_directory+file_name)

    if pname == 'OK_M':
        stnry_df.rename(columns = {'nlat':'lat', 'elon':'lon'}, inplace = True)
    if pname == 'KS_M':
        stnry_df.rename(columns = {'LATITUDE':'lat', 'LONGITUDE':'lon', 'ABBR':'Stn_ID'}, inplace = True)
    if pname == 'ASOS':
        stnry_df.rename(columns = {'CALL':'Stn_ID', 'LAT':'lat', 'LON':'lon'}, inplace = True)

    # Test Data availability
    if d_testing == True:
        try:
            #if there is no files this will will cause the try statement to fail
            stnry_df.iloc[0]
            return True
        except: return False
    # + + + + + + + + + + + + ++ + +

    else: return stnry_df, 'STN_I'

#**************
def read_Radar(config, day, pname, swp=None, rfile= None, d_testing=False, known_scan_time=None):
    #  print("read_Radar(pname=",pname," rfile=", rfile, ")")
    ''' Determine if a given radar is deployed and if so assign the correct location values to it.
    '''
    #  rfile will only be provided if the object being initilized is the main plotting radar
    if rfile != None: 
        ## Det the attributes of the Main Radar (MR)
        if pname == 'WSR88D': # if the main plotting radar is a Wsr88d radar
            index_at_start = rfile.sweep_start_ray_index['data'][swp[0]] #Find the beginning loc of given sweep in data array
            #det the scantime
            MR_time = cftime.num2pydate(rfile.time['data'][index_at_start], rfile.time['units'])
            radar_fields_prep(config, rfile, 'WSR', swp)

        else:# if the main plotting radar is a KA or NOXP radar
            #apply mask and make texture, gradient etc feilds
            if pname == 'NOXP': radar_fields_prep(config, rfile, 'NOXP', swp)
            elif pname in pform_names('KA'): radar_fields_prep(config, rfile, 'KA', swp)
            #det the scantime
            MR_time = datetime.strptime(rfile.time['units'][14:-1], "%Y-%m-%dT%H:%M:%S")

        #det the location info
        MR_lat, MR_lon = rfile.latitude['data'][0], rfile.longitude['data'][0]
        return MR_time, MR_lat, MR_lon , 'MAINR'

    # if no rfile is provided then you are filling in dep info for the Radars not plotting additional data (aka info for marker placement etc)
    elif rfile == None: 
        if pname in pform_names('KA'):
            ## Read in files
            if pname == 'Ka1': r_testing='ka1' # r_testing is the name of the radar you are testing to see if deployed
            elif pname == 'Ka2': r_testing='ka2'
            filename = config.g_TORUS_directory+day+'/data/radar/TTUKa/csv/'+day+'_deployments_'+r_testing+'.csv'

            kadep = []
            #  read in the csv file; if radar didn't dep that day there will be no csv file and the defn will fail (which is ok)
            if os.path.exists(filename): kadep = pd.read_csv(filename)
            else:
                log.debug("    **** File does not exist: %s" %(filename))
                print("File does not exist: %s" %(filename))

            ## If Radar did dep this day det more info about the deployments
            for t in range(kadep.time_begin.count()):
                beginscan = datetime.strptime(kadep.time_begin[t], "%m/%d/%Y %H:%M")
                endscan = datetime.strptime(kadep.time_end[t], "%m/%d/%Y %H:%M") 

                # det if any of the deps occured in the time frame we are interested in: if so record loc and RHI info for the dep
                if known_scan_time >= beginscan and known_scan_time <= endscan:
                    #If defn hasn't failed yet & we entered this if statement then we have dep data relavant to our plot for this radar
                    # If defn was called to det if we had data availability we will exit the defn here
                    if d_testing == True: return True

                    # Otherwise record data about the dep
                    #Record the loc info
                    klat, klon, Head = kadep.lat[t], kadep.lon[t], kadep.heading[t]
                    #Record the RHI info (if applicable... not all dep will have this)
                    try: RHIb, RHIe = kadep.rhib[t], kadep.rhie[t]
                    except: RHIb, RHIe = np.nan, np.nan

                    #set up a namedtuple object to hold the new info
                    r_loc = namedtuple('r_loc', ['lat', 'lon', 'head', 'rhib', 'rhie'])
                    loc = r_loc(lat=klat, lon=klon, head=Head, rhib= RHIb, rhie=RHIe)
                    return loc, 'KA'

        elif pname == 'NOXP':
            # Determine whether the NOXP radar was in the plotting domain (if its not the radar being plotted)
            if rfile == None:
                path = config.g_TORUS_directory+ day+'/data/radar/NOXP/'+day+'/*/sec/*'
                NOXPfiles = sorted(glob.glob(path))
                file = NOXPfiles[0]

                checking=0
                for file in NOXPfiles:
                    head_tail= os.path.split(file)
                    str_time = datetime.strptime(head_tail[1][13:28], "%Y%m%d_%H%M%S")
                    #  str_time = datetime.strptime(head_tail[1][24:39], "%Y%m%d_%H%M%S")
                    end_time = datetime.strptime(head_tail[1][43:58], "%Y%m%d_%H%M%S")
                    valid_time = time_in_range(str_time, end_time, known_scan_time)

                    if valid_time == True:
                        if d_testing == True: return True
                        r = pyart.io.read_cfradial(file)
                        Nlat, Nlon = r.latitude['data'][0], r.longitude['data'][0]

                        #set up a namedtuple object to hold the new info
                        r_loc = namedtuple('r_loc', ['lat', 'lon'])
                        loc = r_loc(lat=Nlat, lon=Nlon)
                        
                        checking= checking+1
                        if checking > 1: print('HEYYYYYYYYY !!!!')

                if checking == 0:
                    if d_testing == True: return False
                return loc, 'NOXP'

        #  Determine what other WSR sites could fall withing the plotting domain
        elif pname == 'WSR88D':
            #  save the nexrad locations to an array from the PyART library
            wsr_locs = pyart.io.nexrad_common.NEXRAD_LOCATIONS
            WSR_df= pd.DataFrame.from_dict(wsr_locs, orient= 'index')
            WSR_df.reset_index(level=None, drop=False, inplace=False, col_level=0, col_fill='')
            WSR_df.index.name = 'R_Name'
            WSR_df.reset_index(inplace= True, drop=False)

            # Test Data availability
            if d_testing == True:
                try: #  if no sites in domain this will cause a failure
                    WSR_df.iloc[0]
                    return True
                except:return False
            else:  return WSR_df, 'WSR'

##########################################################################################
##########
# Classes
##########
class Platform:
    #vars defined in this block (until ######) reamain constant for any object initilized via calling Platform or any Platform subclass
        #  self.var can be retreived latter via typing obj.day etc
        #  vars without self. can be used within Platform or Platform subclasses methods but not for external retrieval
    ######
    def __init__(self, config, day, Name):
        '''setting the ojects attr, __init__ defns will only be run once per object (and will be unique for each object)
        '''
        self.config = config
        self.name = Name #Instance variables
        #get style info for the platforms marker
        self.m_style, self.m_color, self.m_size, self.l_color, self.leg_str, self.leg_entry = pform_attr(self.name)

    # * * *
    @classmethod
    def test_data(self, config, day, pname, scan_time=None):
        '''test if there is a datafile for the platform
        Inputs: pname= the name of the file being tested
        @classmethod allows you to call the defn without creating an object yet via Platforms.test_data
        '''
        try:
            if pname in pform_names('TInsitu'):
                data_avail = read_TInsitu(config, day, pname, d_testing=True)
            elif pname in pform_names('STN_I'):
                data_avail = read_Stationary(config, day, pname, d_testing=True)
            elif pname in pform_names('RADAR'):
                data_avail = read_Radar(config, day, pname, d_testing=True, known_scan_time=scan_time)
            return data_avail
        except:
            error_printing(config.e_test)
            return False

    # * * *
    def getLocation(pform, offsetkm, scan_time=None ,given_bearing= False, pnt_x=None, pnt_y=None):
        ''' This definition has two functions:
                1) If no bearing is specified it will return a namedtuple containing the max/min lat/lons
                    to form a square surrounding the point indicated by lat1,lon1 by x km.
                2) If a bearing is given then the defintion will return one set of lat/lon values
                     (end_lat and end_lon) which is the location x km away from the point lat1, lon1
                     if an observer moved in a straight line the direction the bearing indicated.
        ----
        INPUTS: offsetkm: distance traveled from the starting point to the new locations
                given_bearing: True/False, are you given a specified direction of "travel"
        '''
        def calc(bearing, offsetkm, R, lat1, lon1):
            rad_brng = math.radians(bearing)
            new_lat = np.arcsin(np.sin(lat1) * np.cos(offsetkm/R) + np.cos(lat1) * np.sin(offsetkm/R) * np.cos(rad_brng))
            new_lon = lon1 + np.arctan2(np.sin(rad_brng) * np.sin(offsetkm/R) * np.cos(lat1), np.cos(offsetkm/R) - np.sin(lat1) * np.sin(new_lat))
            return math.degrees(new_lat), math.degrees(new_lon)

        # + + + + + + + + + + + + ++ + +
        #  determine starting point lat and lon
        if pform =='ClickPoint':
            start_lat, start_lon = pnt_y, pnt_x
        
        elif isinstance(pform, Radar):
            start_lat, start_lon = pform.lat, pform.lon

        else:
            # locate the data entry that has the closest time as the plotting radar scantime
            ScanT_entry =pform.df.loc[pform.df['datetime'].sub(scan_time).abs().idxmin()]
            start_lat, start_lon = ScanT_entry.lat, ScanT_entry.lon

            #  if len(start_lat) > 1: start_lat = start_lat.median()
            #  if len(start_lon) > 1: start_lon = start_lon.median()
            if not isinstance(start_lat, np.float32):
                if not isinstance(start_lat, np.float64): start_lat = start_lat.median()
            if not isinstance(start_lon, np.float32): 
                if not isinstance(start_lon, np.float64): start_lon = start_lon.median()

        lat1, lon1 = start_lat * np.pi/180.0 , start_lon * np.pi/180.0
        R = 6378.1 #earth radius (R = ~ 3959 MilesR = 3959)

        # + + + + + + + + + + + + ++ + +
        #  Calc end point(s) lat and lons
        if given_bearing == False:
            for card_dir in [0, 90, 180, 270]:
                new_lat, new_lon = calc(card_dir, offsetkm, R, lat1, lon1)
                
                if card_dir == 0: max_lat= new_lat
                elif card_dir == 90: max_lon= new_lon
                elif card_dir == 180: min_lat= new_lat
                elif card_dir == 270: min_lon= new_lon

            #set up a namedtuple object to hold the new info
            box_extent = namedtuple('box_extent', ['xmin', 'xmax','ymin','ymax'])
            box = box_extent(xmin= min_lon, xmax= max_lon, ymin= min_lat, ymax= max_lat)
            return box

        else: #if a bearing is provided
            end_lat, end_lon = calc(given_bearing, offsetkm, R, lat1, lon1)
            return end_lat, end_lon

    # * * *
    def grab_pform_subset(self, pform, Data, bounding=None, time_offset=None):
        ''' This def will take a given point or pandas dataframe (df) and subset it either spatially or temporially
                1) If time_offset is given the data will be subset temporally
                        Grabs the observed thermo data +/- x mins around radar scan_time.
                        ie if time_offest=5 the resultant subset will span 10 min.
                2) If bounding is given the data will be subset spacially
                        will return the data points that fall within the box defined by ymin,ymax,xmin,xmax
                3) If both scan_time and a bounding region are provided the dataset should be subset both
                        temporally and spatially. (This has not been tested yet so double check if using)
        -----
        INPUTS: Dataframe= True if you are subseting a pandas dataframe, leave as False if you are simply checking wether a point matches the domain
        Returns: The subset dataset (df_sub) and a True/False statement regarding any data in the original dataset
                    matched the subsetting criteria (p_deploy)
        '''
        if pform.type == 'KA' or pform.type == 'NOXP': Single_Point = True
        else: Single_Point = False
        # + + + + + + + + + + + + ++ + +
        #Temporal subset
        ##### + + + + + +
        if time_offset != None:
            if Single_Point == True:  print('Code not written yet to spatially subset a platform without a pandas df (like radars)')
            else:
                lower_tbound = (Data['P_Radar'].Scan_time-dt.timedelta(minutes=time_offset)) 
                upper_tbound = (Data['P_Radar'].Scan_time+dt.timedelta(minutes=time_offset))
                df_sub = pform.df.loc[(pform.df['datetime'] >= lower_tbound) & (pform.df['datetime'] <= upper_tbound)]
                if self.config.print_long == True: print('Dataset has been temporally subset')

        # + + + + + + + + + + + + ++ + +
        #Spatial Subset
        ##### + + + + +
        if bounding != None:
            if Single_Point == True:
                if np.logical_and(pform.lat > bounding.ymin, np.logical_and(pform.lat < bounding.ymax,
                       np.logical_and(pform.lon > bounding.xmin, pform.lon < bounding.xmax))): p_deploy = True
                else: p_deploy = False
            else:
                # if the dataset has not been temporally subset (aka time_offset is none) the dataframe
                #  to do the spatiol subsetting is equivilant to the full original dataset;
                #  otherwise use the previously det df_sub that results from the temporal subset
                if time_offset == None: df_sub = pform.df

                #conduct the spatial subset
                df_sub = df_sub.loc[(df_sub['lat'] >= bounding.ymin) & (df_sub['lat'] <= bounding.ymax) &
                                    (df_sub['lon'] >= bounding.xmin) & (df_sub['lon'] <= bounding.xmax)]
            if self.config.print_long == True: print('Dataset has been spatially subset')

        # + + + + + + + + + + + + ++ + +
        #Determine what to return
        if Single_Point == True: df_sub= 'filler'
        #Test to ensure that there is valid data in the subrange
        #  (aka the pform was dep at the time of Rscan and was within the domain of the plot)
        else:
            try:
                df_sub.iloc[0]
                p_deploy = True
            except:
                p_deploy = False
                error_printing(self.config.e_test)
        return df_sub, p_deploy
    # * * *

####
class Torus_Insitu(Platform):
    def __init__(self, config, day, Name):
        #read in and intilize the dataset; type is the platform type (aka NSSL, UNL, UAS)
        self.df, self.type = read_TInsitu(config, day, Name)
        self.df.name ='{}_df'.format(Name) #assign a name to the pandas dataframe itself
        #add a max/min value, if the object has type NSSL it will apply a mask (this can be easily changed/ mask applied to other platforms)
        if self.type == 'NSSL': mask=True
        else: mask=False
        self.Tv_Min, self.Tv_Max= self.min_max('Thetav', mask)
        self.Te_Min, self.Te_Max= self.min_max('Thetae', mask)
        self.Tf_Min, self.Tf_Max= self.min_max('tfast', mask)
        Platform.__init__(self, config, day, Name)
 
    # * * *
    def min_max(self, var, mask):
        #Mask the dataset(if applicable) and determine the max and min values of pvar:
        #  Masking is based of QC flags etc and can be diff for each platform
        if self.type == "NSSL":
            if mask == True: #mask is a bool, sent to False to by default
                mask_allqc_df = np.ma.masked_where(self.df['all_qc_flags'].values > 0, self.df[var].values) #add masked dataset to the object
                #  self.ts_mask_df = np.ma.masked_where(self.df['qc3'].values > 0, self.df[config.var].values) #add masked dataset to the object
                self.Thetae_ts_mask_df = np.ma.masked_where(self.df['qc3'].values > 0, self.df['Thetae'].values) #add masked dataset to the object
                self.Thetav_ts_mask_df = np.ma.masked_where(self.df['qc3'].values > 0, self.df['Thetav'].values) #add masked dataset to the object
                self.tfast_ts_mask_df = np.ma.masked_where(self.df['qc3'].values > 0, self.df['tfast'].values) #add masked dataset to the object
                if len(mask_allqc_df) == 0: 
                    self.Min, self.Max = np.nan, np.nan
                else:
                    self.Min, self.Max = mask_allqc_df.min(), mask_allqc_df.max()#use the masked dataset to find max and min of var
            elif mask == False:
                self.Min, self.Max = self.df.min()[var], self.df.max()[var] #use the orig dataset to find max/min

        elif self.type == "UNL":
            if mask == True: print('UNL masking code not written yet')
            elif mask == False:
                self.Min, self.Max = self.df.min()[var], self.df.max()[var] #use the orig dataset to find max/min

        elif self.type == "UAS":
            if mask == True: print('UAS masked code not written yet')
            elif mask == False: print('UAS unmasked code not written yet')

        if self.Min== None: self.Min = np.nan
        if self.Max== None: self.Max = np.nan
        return self.Min, self.Max

####
class Stationary_Insitu(Platform):
    def __init__(self, config, day, Name):
        #read in and intilize the dataset;
        self.df, self.type = read_Stationary(config, day, Name)
        self.df.name = '{}_df'.format(Name) #assign a name to the pandas dataframe itself
        Platform.__init__(self, config, day, Name)

####
class Radar(Platform):
    def __init__(self, config, day, Name, Data=None, Rfile= None, Swp=None, Plotting_Radar= False):
        self.config = config
        ## If you are initializing the Plotting Radar object
        if Plotting_Radar == True:
            ## Det the key attr of the main plotting radar and define the class var Scan_time for all objects of Platform
            self.rfile, self.name, self.swp = Rfile, Name, Swp
            self.Scan_time, self.lat, self.lon, self.type= read_Radar(config, day, self.name, self.swp, rfile=self.rfile)
            # Convert time into fancy date string to use in overall plot title
            self.fancy_date_str = self.Scan_time.strftime('%Y-%m-%d %H:%M UTC')

            if self.name in pform_names('KA'):  self.site_name, self.dir_name = self.name, 'KA'
            if self.name == 'WSR88D':  self.site_name, self.dir_name = self.rfile.metadata['instrument_name'], 'WSR'
            if self.name == 'NOXP':  self.site_name, self.dir_name = 'NOXP', 'NOXP'

        ## Otherwise you are initlizing info for radars that doent have to do with plotting actual radar data (aka loc of radars etc but not the data itself)
        elif Plotting_Radar == False:
            Rloc, self.type = read_Radar(config, day, Name, known_scan_time=Data['P_Radar'].Scan_time)

            if self.type == 'KA':
                #store the loc info for a ka radar object
                self.lat, self.lon, self.head, self.rhib, self.rhie = Rloc.lat, Rloc.lon, Rloc.head, Rloc.rhib, Rloc.rhie
            elif self.type == 'NOXP':
                self.lat, self.lon = Rloc.lat, Rloc.lon
            elif self.type == 'WSR':
                #Rloc is a dataframe containing the lats and lons of each WSR88D that is within the plotting domain
                self.df = Rloc
                self.df.name = '{}_df'.format(Name) #assign a name to the pandas dataframe itself
        Platform.__init__(self, config, day, Name)
    
    def plot_topo(self, PLOT, ax):
        if self.config.overlays['GEO']['country_roads'] == True:
            ox.config(log_file=True, log_console=True, use_cache=True) #the config in this line has nothing to do with config.py
            G = ox.graph_from_bbox(PLOT.Domain.ymax, PLOT.Domain.ymin, PLOT.Domain.xmax, PLOT.Domain.xmin)
            ox.save_load.save_graph_shapefile(G, filename='tmp'+str(0), folder=self.config.g_roads_directory , encoding='utf-8')
            fname = self.config.g_roads_directory + 'tmp'+str(0)+'/edges/edges.shp'
            shape_feature = ShapelyFeature(Reader(fname).geometries(), PLOT.trans_Proj, edgecolor='gray', linewidth=0.5)
            ax.add_feature(shape_feature, facecolor='none')
            shutil.rmtree(self.config.g_roads_directory+'tmp'+str(0)+'/')
        if self.config.overlays['GEO']['hwys'] == True:
            fname = self.config.g_roads_directory+'GPhighways.shp'
            shape_feature = ShapelyFeature(Reader(fname).geometries(), PLOT.trans_Proj, edgecolor='grey')#edgecolor='black')
            ax.add_feature(shape_feature, facecolor='none')
        if self.config.overlays['GEO']['county_lines'] == True:
            fname = self.config.g_roads_directory+'cb_2017_us_county_5m.shp'
            shape_feature = ShapelyFeature(Reader(fname).geometries(), PLOT.trans_Proj, edgecolor='gray')
            ax.add_feature(shape_feature, facecolor='none', linewidth=1.5, linestyle="--")
        if self.config.overlays['GEO']['state_lines'] == True:
            states_provinces = cartopy.feature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='10m', facecolor='none')
            ax.add_feature(states_provinces, edgecolor='black', linewidth=2)

