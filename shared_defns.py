#import needed modules
######################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as PathEffects
import pandas as pd
import datetime as dt
from datetime import datetime, date, timedelta
import metpy
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.crs as ccrs
import xarray as xr
from collections import namedtuple
import pyart, sys, traceback, glob, cmocean

## Imports form other files 
############################
import config #this is the file with the plotting controls to access any of the vars in that file use config.var 

#rename a few commonly used vars so that the config.var does not have to be used repeatedly 
print_long, e_test, p_var = config.print_long, config.e_test, config.p_var


################################################################################################
##################
# TROUBLE SHOOTING
##################
def error_printing(e_test):
    ''' Basically I got sick of removing/replacing the try statements while troubleshooting 
    '''
    if e_test == True:
        e = sys.exc_info()
        traceback.print_tb(e[-1])
        #print( "<p>Error: %s</p>" % e )
        print("Error:", e[:-2],'\n')

################################################################################################
###########
# Plotting
###########

################################################################################################
###########
# Data Prep
###########
def pform_names(Type):
    ''' Returns lists of platform names for various platform types; Saves typing in the end 
        ---
        INPUT: Type [str]: valid options include 'ALL','RADAR', 'TInsitu','UNL','NSSL','KA'
        Output: P_list: list of pnames containing the names of the requested type
    '''
    if Type == 'ALL': P_list = ['WTx_M','OK_M','IA_M','ASOS','AWOS','METAR','FFld','LIDR','Prb1','Prb2','WinS','CoMeT1','CoMeT2','CoMeT3','UAS', 'Ka1','Ka2','NOXP','WSR88D']
    elif Type == "STN_I": P_list= ['WTx_M','OK_M','IA_M','ASOS','AWOS', 'METAR']
    elif Type == 'TInsitu': P_list = ['FFld','LIDR','Prb1','Prb2','WinS','CoMeT1','CoMeT2','CoMeT3','UAS']
    elif Type == 'RADAR': P_list = ['Ka1','Ka2','NOXP','WSR88D']
    elif Type == 'MESO': P_list = ['WTx_M','OK_M','IA_M']
    elif Type == 'NWS': P_list= ['ASOS','AWOS','METAR']
    elif Type == 'UNL': P_list = ['CoMeT1','CoMeT2','CoMeT3']
    elif Type == 'NSSL': P_list = ['FFld','LIDR','Prb1','Prb2','WinS']
    elif Type == 'KA': P_list = ['Ka1','Ka2']
    else: print('Please enter a valid list name')
    return P_list

# *******************
def pform_attr(pname, print_long):
    ''' INPUTS: pname: name of the platform
        ----
        RETURNS: legend_entry: The specs regarding the legend presentation for each pform (will be appended latter to make the legend)
                 m_color,m_style etc : the specifications regarding platform presentation on the plots (markershape, color, etc)
    '''
    if print_long== True: print('Made it into pform_attr')

    ##assign the atributes for each platform
    if pname == "Prb1": marker_style, marker_color, line_color, legend_str= '1','xkcd:lightblue','steelblue','Prb1'
    elif pname == "Prb2": marker_style, marker_color, line_color, legend_str= '1','xkcd:watermelon','xkcd:dusty red','Prb2'
    elif pname == "FFld": marker_style, marker_color, line_color, legend_str= '1','xkcd:bubblegum pink','xkcd:pig pink','FFld'
    elif pname == "LIDR": marker_style, marker_color, line_color, legend_str= '1','xkcd:pastel purple','xkcd:light plum','LIDR'
    elif pname == "WinS": marker_style, marker_color, line_color, legend_str= '1','xkcd:peach','xkcd:dark peach','WinS'
    elif pname == "CoMeT1": marker_style, marker_color, line_color, legend_str= '1','brown','brown','CoMeT1'
    elif pname == "CoMeT2": marker_style, marker_color, line_color, legend_str= '1','yellow','yellow','CoMeT2'
    elif pname == "CoMeT3": marker_style, marker_color, line_color, legend_str= '1','black','black','CoMeT3'
    elif pname == "WTx_M": marker_style, marker_color, line_color, legend_str= r'$\AA$' ,'black','black','WTxM Site'   
        #U+1278, #u20a9
    elif pname == 'Ka2': marker_style, marker_color, line_color, legend_str= '8','xkcd:crimson','xkcd:crimson','Ka2'
    elif pname == 'Ka1': marker_style, marker_color, line_color, legend_str= '8','mediumseagreen','mediumseagreen','Ka1'
    elif pname == 'WSR88D': marker_style, marker_color, line_color, legend_str= r'$\Omega$', 'white','black', "WSR88D"

    ##create the legend elements that will be appended to the legend array    
    if pname in pform_names('KA'): #the legend entries for the KA radars
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor='black',markeredgewidth=3,label=legend_str,markerfacecolor=marker_color, markersize=26)
    elif pname == 'WSR88D':
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor=marker_color,markeredgewidth=3,label=legend_str, markersize=26,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')])
    elif pname == 'WTx_M':
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor=marker_color,label=legend_str, markersize=26)
    else:
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor=marker_color,markeredgewidth=3,label=legend_str, markersize=26,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')])
        
    if print_long== True: print('Made it through pform_attr')
    return marker_style,marker_color,line_color,legend_str,legend_entry

#**************
def read_TInsitu(pname, print_long, e_test, tstart=None, tend=None, d_testing=False):
    ''' Reads data files provided on TORUS19 EOL site (Important that datafiles and filenames follow TORUS19 readme)
        ---
        INPUT: filename string (following readme conventions for each platform) 
               tstart,tend: datetime objects
        OUTPUT: dataframe containing all data collected during desired times
    '''
    if pname in pform_names('UNL'): 
        mmfile= glob.glob(config.temploc+config.day+'/mesonets/UNL/UNL.'+pname+'.*')
        mtest=mmfile[0] #if there is no files this will will cause the script to fail (in a good way)
        #if the code has not failed by this point there is data present; if calling defn in a testing capacity 
        #  the defn will exit at this point (saving computing time) otherwise the defn will cont
        if d_testing==True: return True
        # + + + + + + + + + + + + ++ + +
 
        data_hold=[] #empty list to append to 
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
            if tstart is None: hstart = pd_unl['datetime'].iloc[0]
            else: hstart = tstart
            if tend is None: hend= pd_unl['datetime'].iloc[-1]
            else: hend = tend
            # Save only desired iop data
            data_u = pd_unl.loc[(pd_unl['datetime']>=hstart)&(pd_unl['datetime']<=hend)]
            data_hold.append(data_u)
        
        #convert the list holding the dataframes to one large dataframe 
        data_unl=pd.concat(data_hold)
        return data_unl, 'UNL'

    # * * * 
    elif pname in pform_names('NSSL'): 
        mmfile=glob.glob(config.filesys+'TORUS_Data/'+config.day+'/mesonets/NSSL/'+pname+'_'+config.day[2:]+'_QC_met.dat')
        file=mmfile[0]
        #if the code has not failed by this point there is data present; if calling defn in a testing capacity 
        #  the defn will exit at this point (saving computing time) otherwise the defn will cont
        if d_testing==True: return True
        # + + + + + + + + + + + + ++ + +

        # Read NSSL file using column names from readme
        column_names = ['id','time','lat','lon','alt','tfast','tslow','rh','p','dir','spd','qc1','qc2','qc3','qc4']
        data = pd.read_csv(file,header=0,delim_whitespace=True,names=column_names)
        data = data.drop_duplicates()

        # Find timedelta of hours since start of iop (IOP date taken from filename!)
        tiop = dt.datetime(2019, np.int(file[-15:-13]),np.int(file[-13:-11]),0,0,0)
        
        if tstart is None: 
            hstart = tiop
            hstart_dec=hstart.hour+(hstart.minute/60)+(hstart.second/3600) #convert to decimal hours HH.HHH
        else: 
            hstart = (tstart-tiop).seconds/3600
            hstart_dec=hstart

        if tend is None: 
            hend = data['time'].iloc[-1]
        else: 
            hend = (tend-tiop)
            if hend >= dt.timedelta(days=1):  hend = (tend-tiop).seconds/3600 + 24.
            else:  hend = (tend-tiop).seconds/3600
        
        # Save only desired iop data
        data_nssl = data.loc[(data['time']>=hstart_dec)&(data['time']<=hend)]

        # Convert time into timedeltas
        date = dt.datetime.strptime('2019-'+file[-15:-13]+'-'+file[-13:-11], '%Y-%m-%d')
        time_deltas = []
        for i in np.arange(len(data_nssl)):
            j = data_nssl['time'].iloc[i]
            time_deltas = np.append(time_deltas,date+dt.timedelta(hours=int(j),minutes=int((j*60) % 60),seconds=int((j*3600) % 60)))
        data_nssl['datetime']=time_deltas

        ## Caclulate desired variables
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

        t=[]
        for i in range(0, len(data_nssl)):
            t.append([i] *30)
        data_nssl['group'] = tuple(t)
        data_nssl['wmax']=data_nssl.groupby(['group'])['spd'].transform(max)
        return data_nssl, 'NSSL'
    
    # * * * 
    elif pname == 'UAS': 
        if print_long == True: print("no code for UAS yet")
        if d_testing==True: return False
        # + + + + + + + + + + + + ++ + +
	return 'UAS' 

#**************
def read_Stationary(pname, p_domain, print_long, e_test, d_testing=False):
    ''' Determine if a there are any sites from the stationary arrays that fall within the plot domain 
        if so record locations and names of the sites
    '''
    if pname =='WTM':
        meso_df=pd.read_csv(filesys+'radarproc/West_TX_mesonets.csv')   
        data_wtm, pdeploy= grab_pform_subset(print_long, e_test, Data, d_test_df= meso_df, bounding= p_domain, Dataframe= True)
        #if no sites in domain this will cause a failure for d_testing 
        p_test= data_wtm.iloc[1]
        #  if testing for data availability (and the defn has not failed yet) the func will end here
        if d_testing==True: return True
        else: return data_wtm, 'WTM' 

#**************
def radar_dep(pname, scantime, print_long, e_test, d_testing=False): 
    ''' Determine if a given radar is deployed and if so assign the correct location values to it.
    '''
    if pname in pform_names('KA'): 
        ## Read in files 
        if pname == 'Ka1': r_testing='ka1' # r_testing is the name of the radar you are testing to see if deployed
        elif pname == 'Ka2': r_testing='ka2'
        #  read in the csv file; if the Ka radar didn't deploy there will be no csv file and the defn will fail
        kadep=pd.read_csv(config.filesys+'TORUS_Data/'+config.day+'/radar/TTUKa/csv/'+config.day+'_deployments_'+r_testing+'.csv')
        #  if testing for data availability (and the defn has not failed yet) the func will end here
        if d_testing==True: return True
        # + + + + + + + + + + + + ++ + +
        
        # If a Ka Radar had deployed det it's location (lat/lon) and the first and last RHI scan angle of a sweep if applicable.
        define_check=0
        for t in range(kadep.time_begin.count()):
            beginscan=datetime.strptime(kadep.time_begin[t], "%m/%d/%Y %H:%M")
            endscan=datetime.strptime(kadep.time_end[t], "%m/%d/%Y %H:%M")

            if scantime >= beginscan and scantime <= endscan:
                try: klat, klon, head = kadep.lat[t], kadep.lon[t], kadep.heading[t]
                except:
                    klat, klon, head = np.nan, np.nan, np.nan
                    error_printing(e_test)

                try: rhib, rhie = kadep.rhib[t], kadep.rhie[t]
                except:
                    rhib, rhie = np.nan, np.nan
                    error_printing(e_test)
            else: rhib, rhie, klat, klon, head = np.nan, np.nan, np.nan, np.nan, np.nan

            # preseve the actual values so they are not overwritten by other deployments
            if str(klat) != 'nan':
                KLAT, KLON, HEAD, RHIB, RHIE = klat, klon, head, rhib, rhie
                define_check = 1

        #In the instance that the second radar isnt deployed at time of radarscan the KLAT, KLON,.... will still be defined
          ## If the other radar is deployed it will not overwrite the assigned values
        if define_check == 0: KLAT, KLON, HEAD, RHIB, RHIE = np.nan, np.nan, np.nan, np.nan, np.nan

        #set up a namedtuple object to hold the new info
        r_loc = namedtuple('r_loc', ['lat','lon','head','rhib','rhie'])
        loc =r_loc(lat=KLAT, lon=KLON, head=HEAD, rhib= RHIB, rhie=RHIE)
        return loc, 'KA' 

    # * * * 
    elif pname == 'WSR88D':
        #save the nexrad locations to an array from the PyART library
        locs = pyart.io.nexrad_common.NEXRAD_LOCATIONS
        print(locs)
        data_WSR88D, pdeploy= grab_pform_subset(print_long, e_test, Data, d_test_df= locs, bounding= p_domain, Dataframe=True)
        #if no sites in domain this will cause a failure for d_testing 
        p_test= data_WSR88D.iloc[1]
        #  if testing for data availability (and the defn has not failed yet) the func will end here
        if d_testing==True: return True
        else: return data_wtm, 'WSR' 

# * * * * *
def Add_to_DATA(DType,Data,subset_pnames,print_long,rfile=None,rname=None,swp=None):
    if print_long== True: print('Made it into Add_to_DATA')
    
    if DType == 'TInsitu':
        #for these pforms mainly interested if we have any valid data for the day not worried about spatial/temporal subset yet  
        #  To save computing time this data will only be read in once and anyfurther subsetting will be done later
        for pname in pform_names(DType):
            data_avail=Platform.test_data(pname) #for each platform test to see if we have data for that day 

            if data_avail==True:
                subset_pnames.append(pname) #append the pname to the subset_pnames list 
                #  load data for the pform (aka initialize an object of the appropriate class); place in dict with key of pname
                Data.update({pname: Torus_Insitu(pname)}) 
                if print_long==True: print("Data can be read in for platform %s" %(pname))
            else: 
                if print_long==True: print("No data available to be read in for platform %s" %(pname))

        #go through the values stored in Data. If any are an object of the Torus_Insitu Class then it will add a max/min value, 
        #  if the object has type NSSL it will apply a mask (this can be easily changed/ mask applied to other platforms)
        for p in Data.values():
            if isinstance(p,Torus_Insitu):
                if p.type == 'NSSL': p.min_max(config.p_var,mask=True)
                else: p.min_max(config.p_var)

    # * * * 
    elif DType == 'STN_I':
        #What sites are contained within our plotting domain can change for each image so will be updated for every plot  
        for pname in pform_names('MESO'):
            if pname == 'WTx_M':
                data_avail=Platform.test_data(pname, p_domain)
                
                if data_avail==True:
                    #  dont repeatedly append to the list if the platform is already included
                    if pname in subset_pnames: pass 
                    else: subset_pnames.append(pname) #append the pname to the subset_pnames list 
                    #  load loc data for the sites in plotting domain
                    Data.update({pname: Stationary_Insitu(pname)}) 
                    if print_long==True: print("Data can be read in for platform %s" %(pname))
                else: 
                    if print_long==True: print("No data available to be read in for platform %s" %(pname))

            elif pname in ['OK_M','IA_M']: print('code not written yet')
   
        if pname in ['ASOS', 'AWOS', 'METAR']: print('code not written yet')

    # * * * 
    elif DType == 'RADAR':
        #this can change for each image (aka from diff scans of the plotting radar) so will be updated for every plot  
        #  ie. are the radars deployed at a given time which (if any) WSR88D fall within plotting domain etc 
        for pname in pform_names('RADAR'):
            if pname in pform_names('KA'):
                data_avail=Platform.test_data(pname, Data)

                if data_avail==True:
                    #  dont repeatedly append to the list if the platform is already included
                    if pname in subset_pnames: pass 
                    else: subset_pnames.append(pname) #append the pname to the subset_pnames list 
                    #  load loc data (and in this case rhi angles if applicable)
                    Data.update({pname: Radar(pname,Data)}) 
                    if print_long==True: print("Data can be read in for platform %s" %(pname))
                else: 
                    if print_long==True: print("No data available to be read in for platform %s" %(pname))

            if pname == 'WSR88D':
                data_avail=Platform.test_data(pname, Data)

                if data_avail==True:
                    #  dont repeatedly append to the list if the platform is already included
                    if pname in subset_pnames: pass 
                    else: subset_pnames.append(pname) #append the pname to the subset_pnames list 
                    #  load loc data (and in this case rhi angles if applicable)
                    Data.update({pname: Radar(pname,Data)}) 
                    if print_long==True: print("Data can be read in for platform %s" %(pname))
                else: 
                    if print_long==True: print("No data available to be read in for platform %s" %(pname))
        
            if pname == 'NOXP': print('code not written yet')
    
    # * * * 
    elif DType == 'Plotting_Radar':
        #this oject contains the data relevant to the particular radar data being plotted for each image 
        #  (aka platform name, scantime, the data file etc) this will be updated for each image
        Data.update({'P_Radar':Radar(rname, Rfile=rfile,Swp=swp,Plotting_Radar=True)})

    # * * * 
    elif DType == 'PVar':
        #Create a Pvar object add it to the data dict and find the global max and min
        #  This object will only be once updated once (aka same for all images for a given run)
        Data.update({'Var': Pvar(config.p_var)})
        Data['Var'].find_global_max_min(Data)
        Data['Var'].make_dummy_plot()

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
    #  getattr(abject,attribute): get the attibute value of the object, if object doesnt have the attr can return a default value
    #  hasattr(object, attribute): like getattr but returns a true/false bool
    #  vars(object): returns all the values of the attributes for a given object
    #  isinstance
    #  issubclass
    #  there are more: check google for built in class functionalitly
    if print_long== True: print('Made it through Add_to_DATA')
    return Data,subset_pnames 


################################################################################################
##########
# Classes
##########
class Platform:
    #vars defined in this block (until ######) reamain constant for any object initilized via calling Platform or any Platform subclass
        #  self.var can be retreived latter via typing obj.day etc 
        #  vars without self. can be used within Platform or Platform subclasses methods but not for external retrieval
    Day= config.day #Class variables
    Tstart, Tend = config.tstart, config.tend
    Print_long, E_test = config.print_long, config.e_test
    P_domain= 'Place holder' #named tuple defining the bounds of the plot domain 
    Scan_time= 'Place holder'
    ######
    def __init__(self, Name):
        '''setting the ojects attr, __init__ defns will only be run once per object (and will be unique for each object)
        '''
        self.name=Name #Instance variables
        #get style info for the platforms marker
        self.m_style, self.m_color, self.l_color, self.leg_str, self.leg_entry= pform_attr(self.name, self.Print_long)
    

    @classmethod
    def test_data(self, pname, Data=None):
        '''test if there is a datafile for the platform
        Inputs: pname= the name of the file being tested
                Data= the dict containing previous ojects (only needed if testing a radar platform)
        @classmethod allows you to call the defn without creating an object yet via Platforms.test_data 
        '''
        try:
            if pname in pform_names('TInsitu'):
                data_avail=read_TInsitu(pname, self.Print_long, self.E_test, self.Tstart, self.Tend, d_testing=True)
            elif pname in pform_names('STN_I'):
                data_avail=read_Stationary(pname, p_domain, self.Print_long, self.E_test, d_testing=True)
            elif pname in pform_names('RADAR'):
                data_avail=radar_dep(pname, Data['P_Radar'].time, self.Print_long, self.E_test, d_testing=True)
            return data_avail
        except: 
            error_printing(e_test)
            return False
        

    def getLocation(self, offsetkm, given_bearing= False):
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
        #  determine starting point lat and lon
        start_lat, start_lon = self.lat, self.lon  
        lat1, lon1 = start_lat * np.pi/180.0 , start_lon * np.pi/180.0
        R = 6378.1 #earth radius (R = ~ 3959 MilesR = 3959)

        if given_bearing == False:
            for brng in [0, 90, 180, 270]:
                bearing = (brng/90.) * np.pi/2.

                new_lat = np.arcsin(np.sin(lat1) * np.cos(offsetkm/R) + np.cos(lat1) * np.sin(offsetkm/R) * np.cos(bearing))
                new_lon = lon1 + np.arctan2(np.sin(bearing) * np.sin(offsetkm/R) * np.cos(lat1), np.cos(offsetkm/R) - np.sin(lat1) * np.sin(new_lat))
                new_lat, new_lon = 180.0 * new_lat/np.pi , 180.0 * new_lon/np.pi

                if brng == 0: max_lat= new_lat
                elif brng == 90: max_lon= new_lon
                elif brng == 180: min_lat= new_lat
                elif brng == 270: min_lon= new_lon

            #set up a namedtuple object to hold the new info
            box_extent = namedtuple('box_extent', ['ymin','ymax','xmin','xmax'])
            box = box_extent(ymin= min_lat, ymax= max_lat, xmin= min_lon, xmax= max_lon)
            return box

        else: #if a bearing is provided
            bearing = (given_bearing/90.) * np.pi/2.

            new_lat = np.arcsin(np.sin(lat1) * np.cos(offsetkm/R) + np.cos(lat1) * np.sin(offsetkm/R) * np.cos(bearing))
            new_lon = lon1 + np.arctan2(np.sin(bearing) * np.sin(offsetkm/R) * np.cos(lat1), np.cos(offsetkm/R) - np.sin(lat1) * np.sin(new_lat))
            end_lat, end_lon= 180.0 * new_lat/np.pi , 180.0 * new_lon/np.pi
            return end_lat, end_lon


    def grab_pform_subset(self, print_long, e_test, Data, bounding=None, time_offset=None, d_test_df=None, Dataframe=False):
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
        #if doing datatest on stationary
        if d_test_df != None:
            self.df=d_test_df

        #Temporal subset
        if time_offset != None:
            if Dataframe== False:  print('code not written yet')
            elif Dataframe== True:
                aaa = self.df.loc[(self.df['datetime'] >= Data['P_Radar'].time-dt.timedelta(minutes=time_offset))]
                df_sub = aaa.loc[(aaa['datetime'] <= Data['P_Radar'].time+dt.timedelta(minutes=time_offset))]
            if print_long == True: print('Dataset has been temporally subset')

        #Spatial Subset
        if bounding != None:
            if Dataframe== False: 
                print(bounding)
                print(type(bounding))
                if np.logical_and(self.lat > bounding.ymin, np.logical_and(self.lat < bounding.ymax, 
                       np.logical_and(self.lon > bounding.xmin, self.lon < bounding.xmax))):
                    p_deploy = True
                else: p_deploy = False
            elif Dataframe== True:
                #if both time_offset and bounding area is given then the spatial subset start from the already
                    #temporally subset dataframe... if not will start with the full platform dataframe
                if time_offset != None: aaa= df_sub.loc[(self.df['lat']>=bounding.ymin)]
                else: aaa = self.df.loc[(self.df['lat']>=bounding.ymin)]

                bbb = aaa.loc[(aaa['lat']<=bounding.ymax)]
                ccc = bbb.loc[(bbb['lon']>=bounding.xmin)]
                df_sub = ccc.loc[(ccc['lon']<=bounding.xmax)]
            if print_long == True: print('Dataset has been spatially subset')

        #Determine what to return 
        if Dataframe== False: return p_deploy
        elif Dataframe == True:
            #Test to ensure that there is valid data in the subrange
            #  (aka the platform was deployed during the time of radarscan or that any of the points fall within the map area)
            try:
                p_test= df_sub.iloc[1]
                p_deploy = True
            except:
                p_deploy = False
                error_printing(e_test)
            return df_sub, p_deploy

####
class Torus_Insitu(Platform):
    def __init__(self, Name):
        #read in and intilize the dataset; type is the platform type (aka NSSL, UNL, UAS) 
        self.df, self.type= read_TInsitu(Name, self.Print_long, self.E_test, self.Tstart, self.Tend)
        self.df.name='{}_df'.format(Name) #assign a name to the pandas dataframe itself
        Platform.__init__(self, Name)


    def min_max(self,p_var,mask=False):
        #Mask the dataset(if applicable) and determine the max and min values of pvar: 
        #  Masking is based of QC flags etc and can be diff for each platform
        if self.type == "NSSL":
            if mask== True: #mask is a bool, sent to False to by default
                self.mask_allqc_df= np.ma.masked_where(self.df['qc_flag'].values>0, self.df[p_var].values) #add masked dataset to the object
                ts_mask_df = np.ma.masked_where(self.df['qc3'].values>0, self.df[config.p_var].values) #add masked dataset to the object
                self.Min, self.Max = self.mask_allqc_df.min(), self.mask_allqc_df.max()#use the masked dataset to find max and min of p_var
            elif mask == False:
                self.Min, self.Max = self.df.min()[p_var], self.df.max()[p_var] #use the orig dataset to find max/min

        elif self.type == "UNL":
            if mask== True:  print('no code written for this yet')
            elif mask == False:
                self.Min, self.Max = self.df.min()[p_var], self.df.max()[p_var] #use the orig dataset to find max/min

        elif self.type == "UAS":
            if mask== True:  print('no code written for this yet')
            elif mask == False: print('code not written yet')

        else: print("what platform is this ? ", self.type)


    def plot_Tpform(self, Data, ax, print_long, e_test, border_c='xkcd:light grey',labelbias=(0,0)):
        ''' Plot the in situ platform markers, barbs and pathline
        ----
        INPUTS: file: the in situ pandas dataframe
                var: dictionary containing info relating to p_var (ie name, max, min)
                p_attr: dictionary containing platform styling info (color, shape etc)
                ax: axes of the subplot to plot on

        Optional Inputs: border_c, labelbias: color of background for the pathline, if you want to add labels directly to the plot this can offset it from the point
        '''
        if print_long== True: print('made it into platform_plot')

        #grab the subset of data of +- interval around radar scan
        p_sub, p_deploy = self.grab_pform_subset(print_long, e_test, Data, time_offset=config.cline_extent, Dataframe= True)

        if p_deploy == False:
            if print_long==True: print('The platform was not deployed at this time')

        #if there is data for the platform that falls within the time and location of interest
        elif p_deploy == True:
            ##Plot the line that extends +/- min from the platform location; The colorfill indicates values of the specifiec p_var (ie Thetae etc)
            #  fill in the values of the colorfill
            C = cmocean.cm.curl((p_sub[config.p_var].values - Data['Var'].global_min) / (Data['Var'].global_max - Data['Var'].global_min))
            ax = plt.gca()
            for i in np.arange(len(p_sub['lon'])-1): 
                #This is the border of the colorline
                x=p_sub['lon'].values
                y=p_sub['lat'].values
                ax.plot([x[i], x[i+1]], [y[i], y[i+1]], c=border_c, linewidth=10.5, transform=ccrs.PlateCarree(), zorder=3)
                #This is the colorramp colorline
                ax.plot([x[i], x[i+1]], [y[i], y[i+1]], c=C[i], linewidth=7.5, transform=ccrs.PlateCarree(), zorder=4)

            #find the value of the index that is halfway through the dataset (this will be the index associated with radar_scantime)
            mid_point=(p_sub.index[-1]-p_sub.index[0])/2
            col_lon, col_lat =p_sub.columns.get_loc('lon'), p_sub.columns.get_loc('lat')
            col_U, col_V =p_sub.columns.get_loc('U'), p_sub.columns.get_loc('V')

            #plot the platform marker at the time closest to the scantime (aka the time at the halfway point of the subset platform dataframe)
            ax.plot(p_sub.iloc[mid_point,col_lon], p_sub.iloc[mid_point,col_lat], transform=ccrs.PlateCarree(), marker=self.m_style, markersize=20, markeredgewidth='3',color=self.m_color, path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')], zorder=10)
            #plot labels for the marker on the plot itself
            #  d2=plt.text(p_sub.iloc[mid_point,col_lon]+labelbias[0], p_sub.iloc[mid_point,col_lat]+labelbias[1], p.name, transform=ccrs.PlateCarree(), fontsize=20, zorder=9, path_effects=[patheffects.withstroke(linewidth=4,foreground=color)])

            #plot a dot at the end of the colorline in the direction the platform is moving (aka the last time in the subset dataframe)
            ax.plot(p_sub.iloc[-1,col_lon], p_sub.iloc[-1,col_lat],transform=ccrs.PlateCarree(), marker='.',markersize=10, markeredgewidth='3',color='k',zorder=9)

            #plot windbarbs
            #  p_sub.iloc[::x,col_index] returns every x'th value
            stationplot = metpy.plots.StationPlot(ax, p_sub.iloc[::30,col_lon], p_sub.iloc[::30,col_lat], clip_on=True, transform=ccrs.PlateCarree())
            stationplot.plot_barb(p_sub.iloc[::30,col_U], p_sub.iloc[::30,col_V],sizes=dict(emptybarb= 0),length=7)
        if print_long== True: print('made it through platform_plot')

####
class Stationary_Insitu(Platform):
    def __init__(self, Name):
        #read in and intilize the dataset; 
        self.df, self.type= read_Stationary(Name, p_domain, self.Print_long, self.E_test)
        self.df.name='{}_df'.format(Name) #assign a name to the pandas dataframe itself
        Platform.__init__(self, Name)
                    
####
class Radar(Platform):
    #Class variable: pyart display object for the Main plotting radarfile (will be defined after MainRADAR is det)  
    display= 'Place Holder'

    def __init__(self, Name, Data=None, Rfile= None, Swp=None, Plotting_Radar= False):
        if Plotting_Radar == False:
            Rloc, self.type = radar_dep(Name, Data['P_Radar'].time, self.Print_long, self.E_test)
            self.lat, self.lon= Rloc.lat, Rloc.lon
            if self.type == 'KA':
                self.head, self.rhib, self.rhie=Rloc.head, Rloc.rhib, Rloc.rhie
            elif self.type =='WSR':
                self.WSR_df=Rloc
            Platform.__init__(self, Name)

        if Plotting_Radar==True:
            # Determine the key attributes of the main plotting radar (such as scantime etc)
            self.rfile, self.name, self.swp= Rfile, Name, Swp 
            if self.name in pform_names('KA'):
                self.azimuth = self.rfile.fixed_angle['data'][self.swp]
                self.time= datetime.strptime(self.rfile.time['units'][14:-1], "%Y-%m-%dT%H:%M:%S") #Time at which scanning began
            elif self.name == 'WSR88D':
                index_at_start = self.rfile.sweep_start_ray_index['data'][self.swp] #Find the beginning loc of given sweep in data array
                self.time= num2date(self.rfile.time['data'][index_at_start], self.rfile.time['units']) #Time at which scanning began
            elif self.name == 'NOXP': print('code not written yet')
        
            # Convert time into fancy date string to use in title
            self.fancy_date_str = self.time.strftime('%Y-%m-%d %H:%M UTC')
            self.type='MainRADAR'


    def rhi_spokes_rings(self):
        ''' Plot the RHI spoke and ring for a radar
        '''
        if print_long== True: print('made it into rhi_spokes_rings')
        #produce spoke and ring
        for j in range(int(self.rhib), int(self.rhie)+1, 10):
            ang= self.head + j
            if ang > 360.: ang=int(ang-360.)
            
            #  radius = Data['P_Radar'].rfile.range['data'][-1]-500.)/1000.
            if self.type == 'KA': radius=20.905 
            else: print('code not written for other radars yet')
            
            #this plots a circle that connects the spokes
            latArray,lonArray = [], []
            for bearing in range(int(self.head + self.rhib), int(self.head + self.rhie+1)): #degrees of sector
                lat2, lon2 = getLocation(radius, given_bearing=bearing)
                latArray.append(lat2)
                lonArray.append(lon2)
            self.display.plot_line_geo(lonArray, latArray, marker=None, color='grey', linewidth=.25) #this plots a circle that connects the spokes
           
            #plt the spokes 
            C, D= getLocation(radius, given_bearing=ang)
            self.display.plot_line_geo([self.lon, D], [self.lat, C], marker=None, color='k', linewidth=0.5, linestyle=":")

            ## optional labels
            #if np.logical_and(C>ymin,np.logical_and(C<ymax,np.logical_and(D>xmin,D<xmax))):
                    #d1=plt.text(D, C, str(ang),horizontalalignment='center',transform=ccrs.PlateCarree(),fontsize=10,zorder=9,path_effects=([PathEffects.withStroke(linewidth=4,foreground='xkcd:pale blue')])
        if print_long== True: print('made it through rhi_spokes_rings')

#########################################
### set up Pvar class (this is not a subclass of Platform)
class Pvar:
    def __init__(self,p_var):
        self.name=p_var
        #establish lable for the colorbar and tseries ylabel 
        if self.name == "Thetae": self.v_lab= "Equivalent Potential Temp [K]"
        elif self.name == "Thetav": self.v_lab= "Virtual Potential Temp [K]"

    def find_global_max_min(self, Dict):
        '''determine the global max and min across all platforms for p_var
        '''
        val_hold = []
        for p in Dict.values():
            if hasattr(p,'Min')==True: val_hold.append(p.Min)
            if hasattr(p,'Max')== True:val_hold.append(p.Max)
        self.global_min, self.global_max=min(val_hold), max(val_hold)

    def make_dummy_plot(self):
        ''' Make a dummy plot to allow for ploting of colorbar
        '''
        cmap=cmocean.cm.curl
        Z = [[0,0],[0,0]]
        levels = np.arange(self.global_min,self.global_max+1,1)
        CS3 = plt.contourf(Z, levels, cmap=cmap)#,transform=datacrs)
        plt.clf()
        self.CS3= CS3
######################################
        
### set up PLOT class (this is not a subclass of Platform or Pvar)
#  class PLOT:
    #  @staticmethod
    #  def plot_bground_features():
