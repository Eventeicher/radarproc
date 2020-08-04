import numpy as np
from matplotlib.lines import Line2D
import matplotlib.patheffects as PathEffects
import pandas as pd
import datetime as dt
from datetime import timedelta
import metpy.calc as mpcalc
from metpy.units import units
import xarray as xr
import sys, traceback
import glob
import warnings
warnings.filterwarnings(action='ignore')

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
################################################################################################

#**************
def read_platforms(pname,day,print_long,e_test,tstart=None,tend=None,d_testing=False):
    ''' Reads data files provided on TORUS19 EOL site
        Important that datafiles and filenames follow TORUS19 readme
        INPUT: filename string (following readme conventions for each platform) 
               tstart,tend: datetime objects
        OUTPUT: dataframe containing all data collected during desired times
    '''
    if pname in ['CoMeT1', 'CoMeT2', 'CoMeT3']:
        mmfile= glob.glob('/Volumes/Samsung_T5/Research/TORUS_Data/'+day+'/mesonets/UNL/UNL.'+pname+'.*')
        #if there is no files this will will cause the script to fail (in a good way)
        mtest=mmfile[0]
        #if the code has not failed by this point there is data present; if calling defn in a testing capacity the defn will exit at this point (saving computing time) otherwise the defn will cont
        if d_testing==True:
                return True
 
        #empty list to append to 
        data_hold=[]
        for i in range(len(mmfile)):
            ds = xr.open_dataset(mmfile[i])
            
            #convert form epoch time to utc datetime object
            timearray = np.array([dt.datetime.utcfromtimestamp(t/1e9) for t in ds.time.values])
            U,V = mpcalc.wind_components(ds.wind_speed, ds.wind_dir)

            lats = np.array([40])*units.degrees
            #create new xarray dataset to make plotting easier cause original is trash
            dims = ['datetime']
            coords = {
                'datetime': timearray
            }
            data_vars = {
                'lat': (dims, ds.lat.values, {'units':str(lats.units)}),
                'lon': (dims, ds.lon.values, {'units':str(lats.units)}),
                'Z_ASL': (dims, ds.alt.values, {'units':str(ds.alt.units)}),
                'Z_AGL': (dims, np.zeros_like(ds.alt.values), {'units':str(ds.alt.units)}),
                'Temperature': (dims, ds.fast_temp.values, {'units':str(ds.fast_temp.units)}),
                'Dewpoint': (dims, ds.dewpoint.values, {'units':str(ds.dewpoint.units)}),
                'RH': (dims, ds.calc_corr_RH.values, {'units':str(ds.calc_corr_RH.units)}),
                'Pressure': (dims, ds.pressure.values, {'units':str(ds.pressure.units)}),
                'U': (dims, U.m, {'units':str(U.units)}),
                'V': (dims, V.m, {'units':str(V.units)}),
                'Theta': (dims, ds.theta.values, {'units':str(ds.theta.units)}),
                'Thetav': (dims, ds.theta_v.values, {'units':str(ds.theta_v.units)}),
                'Thetae': (dims, ds.theta_e.values, {'units':str(ds.theta_e.units)})
            }
            subds = xr.Dataset(data_vars, coords)
            
            #convert to pandas 
            pd_unl=subds.to_dataframe()
            pd_unl.reset_index(inplace=True)

            # Subset the dataset to the desired time
            if tstart is None:
                hstart = pd_unl['datetime'].iloc[0]
            else:
                hstart = tstart

            if tend is None:
                hend= pd_unl['datetime'].iloc[-1]
            else:
                hend   = tend
            
            # Save only desired iop data
            data_u = pd_unl.loc[(pd_unl['datetime']>=hstart)&(pd_unl['datetime']<=hend)]
            data_hold.append(data_u)
        
        #convert the list holding the dataframes to one large dataframe 
        data_unl=pd.concat(data_hold)
        return data_unl
    
    elif pname in ['FFld','LIDR','Prb1','Prb2','WinS']:

        mmfile=glob.glob('/Users/severe2/Research/TORUS_data/'+day+'/mesonets/NSSL/'+pname+'_'+day[2:]+'_QC_met.dat')
        file=mmfile[0]
        
        #if the code has not failed by this point there is data present; if calling defn in a testing capacity the defn will exit at this point (saving computing time) otherwise the defn will cont
        if d_testing==True:
            return True
        
        column_names = ['id','time','lat','lon','alt',
                        'tfast','tslow','rh','p','dir',
                        'spd','qc1','qc2','qc3','qc4']
        # Read NSSL file using column names from readme
        data = pd.read_csv(file,header=0,delim_whitespace=True,names=column_names)
        data = data.drop_duplicates()

        # Find timedelta of hours since start of iop
        # (IOP date taken from filename!)
        tiop = dt.datetime(2019, np.int(file[-15:-13]),np.int(file[-13:-11]),0,0,0)
        if tstart is None:
            hstart = tiop
            #convert to decimal hours HH.HHH
            hstart_dec=hstart.hour+(hstart.minute/60)+(hstart.second/3600)
        else:
            hstart = (tstart-tiop).seconds/3600
            hstart_dec=hstart

        if tend is None:
            hend = data['time'].iloc[-1]
        else:
            hend   = (tend-tiop)
            if hend >= dt.timedelta(days=1):
                hend = (tend-tiop).seconds/3600 + 24.
            else:
                hend = (tend-tiop).seconds/3600
        # Save only desired iop data
        data_nssl = data.loc[(data['time']>=hstart_dec)&(data['time']<=hend)]

        # Convert time into timedeltas
        date = dt.datetime.strptime('2019-'+file[-15:-13]+'-'+file[-13:-11], '%Y-%m-%d')
        time_deltas = []
        for i in np.arange(len(data_nssl)):
            j = data_nssl['time'].iloc[i]
            time_deltas = np.append(time_deltas,date+dt.timedelta(hours=int(j),minutes=int((j*60) % 60),seconds=int((j*3600) % 60)))
        data_nssl['datetime']=time_deltas

        # Caclulate desired variables
        p,t   = data_nssl['p'].values*units.hectopascal, data_nssl['tfast'].values*units.degC
        theta = mpcalc.potential_temperature(p,t)
        data_nssl['Theta'] = theta.magnitude

        r_h     = data_nssl['rh'].values/100
        mixing = mpcalc.mixing_ratio_from_relative_humidity(r_h,t,p)
        thetav = mpcalc.virtual_potential_temperature(p,t,mixing)
        data_nssl['Thetav'] = thetav.magnitude

        td     = mpcalc.dewpoint_rh(temperature=t, rh=r_h)
        thetae = mpcalc.equivalent_potential_temperature(p,t,td)
        data_nssl['Thetae'] = thetae.magnitude

        Spd =data_nssl['spd'].values*units('m/s')
        dire =data_nssl['dir'].values*units('degrees')
        u,v =mpcalc.wind_components(Spd,dire)
        data_nssl['U']=u.to('knot')
        data_nssl['V']=v.to('knot')

        q_list=['qc1','qc2','qc3','qc4']
        data_nssl['qc_flag']=data_nssl[q_list].sum(axis=1)

        return data_nssl
    
    elif pname in ['UAS']:
        if print_long == True: print("no code for UAS yet")
        if d_testing==True:
            return False
    elif pname in  ['Ka1', 'Ka2']:
        #read in file
        if pname == 'Ka1':
            ka_files =sorted(glob.glob(filesys+'TORUS_Data/'+day+'/radar/TTUKa/netcdf/ka1/dealiased_*'))
        elif pname == 'Ka2':
            ka_files =sorted(glob.glob(filesys+'TORUS_Data/'+day+'/radar/TTUKa/netcdf/ka2/dealiased_*'))
        #test to see if there is data if not it will fail 
        mtest=ka_file[0]
        
        if d_testing==True:
            return True
        
        for thefile in ka_files[:]:
            if print_long == True: print(str(thefile))
            radar = pyart.io.read(thefile)

        return radar 

# *******************

def platform_attr(pname, print_long):
    ''' Assign attributes such as color, markershape, label etc to each platform 
    ----
    INPUTS:
    p_file: pandas dataframe for the platform
    l_array: Array, each platforms information is appended to this array which is used to plot the legend
                    You should define this variable if you are calling the defn while in radar subplots otherwise leave blank.
                    This prevents a platform being added to the legend twice
    radar_m: True/False, UNL and NSSL require file.name to access their str identifier otherplatforms do not
    rad_site: How I am handeling the labeling of NEXRAD at the moment ..... will prob come back and remove
    r_s : ***********fill in comment 
    ----
    RETURNS:
    legend_elements: Array 
    P_Attr: dict, contains the attibute info for the given platform 
    '''
    if print_long== True: print('Made it into platform_attr')

    ##assign the atributes for each platform
    if pname == "Prb1":
        marker_style, marker_color, line_color, legend_str= '1','xkcd:lightblue','steelblue','Prb1'
    elif pname == "Prb2":
        marker_style, marker_color, line_color, legend_str= '1','xkcd:watermelon','xkcd:dusty red','Prb2'
    elif pname == "FFld":
        marker_style, marker_color, line_color, legend_str= '1','xkcd:bubblegum pink','xkcd:pig pink','FFld'
    elif pname == "LIDR":
        marker_style, marker_color, line_color, legend_str= '1','xkcd:pastel purple','xkcd:light plum','LIDR'
    elif pname == "WinS":
        marker_style, marker_color, line_color, legend_str= '1','xkcd:peach','xkcd:dark peach','WinS'
    elif pname == "CoMeT1":
        marker_style, marker_color, line_color, legend_str= '1','brown','brown','CoMeT1'
    elif pname == "CoMeT2":
        marker_style, marker_color, line_color, legend_str= '1','yellow','yellow','CoMeT2'
    elif pname == "CoMeT3":
        marker_style, marker_color, line_color, legend_str= '1','black','black','CoMeT3'
    elif pname == "WTx Mesonet":
        marker_style, marker_color, line_color, legend_str= r'$\AA$' ,'black','black','WTM Site'   
        #U+1278
        #u20a9
    elif pname == 'KA2':
        marker_style, marker_color, line_color, legend_str= '8','xkcd:crimson','xkcd:crimson','Ka2'
    elif pname == 'KA1':
        marker_style, marker_color, line_color, legend_str= '8','mediumseagreen','mediumseagreen','Ka1'
    elif pname == 'WSR88D':
        marker_style, marker_color, line_color, legend_str= r'$\Omega$', 'white','black', "WSR88D"

    ##create the legend elements that will be appended to the legend array    
    if pname in ['KA1','KA2']: #the legend entries for the KA radars
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor='black',markeredgewidth=3,label=legend_str,markerfacecolor=marker_color, markersize=26)
    elif pname == 'WSR88D':
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor=marker_color,markeredgewidth=3,label=legend_str, markersize=26,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')])
    elif pname == 'WTxM':
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor=marker_color,label=legend_str, markersize=26)
    else:
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor=marker_color,markeredgewidth=3,label=legend_str, markersize=26,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')])
        
    if print_long== True: print('Made it through platform_attr')
    return marker_style,marker_color,line_color,legend_str,legend_entry

# *******************

def maskdata(p_var, platform_file, mask=True):
    ''' Read in a dataset and return the masked version of it 
            Masking is based of QC flags etc and can be diff for each platform 
            If mask= False then the defn will simply return the dataset unmodified
    '''
    platform_unmasked= platform_file[p_var].values
    
    if mask== False:
        platform_data= platform_unmasked
    elif mask== True:
        platform_name= str(platform_file.name)
        if (platform_name in ['FFld_df','WinS_df','LIDR_df','Prb1_df','Prb2_df']):
            platform_data= np.ma.masked_where(platform_file['qc_flag'].values>0, platform_unmasked)
        elif (platform_name in ['CoMeT1_df','CoMeT2_df','CoMeT3_df']):
            #for now
            platform_data= platform_unmasked
        elif (platform_name in ['Insert UAS filenames here']):
            print("will be filled in")
        else:
            print("What platform are you trying to use?")
    
    return platform_data
