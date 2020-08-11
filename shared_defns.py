#import needed modules
######################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as PathEffects
import pandas as pd
import datetime as dt
from datetime import datetime, date, timedelta
import metpy.calc as mpcalc
from metpy.units import units
import xarray as xr
from collections import namedtuple
import sys, traceback, glob, cmocean

## Imports form other files 
############################
import config #this is the file with the plotting controls to access any of the vars in that file use config.var 

#rename a few commonly used vars so that the config.var does not have to be used repeatedly 
print_long=config.print_long
e_test=config.e_test
p_var=config.p_var
filesys=config.filesys
temploc=config.temploc
            

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

    else: #if a bearing is provided
        bearing = (given_bearing/ 90.)* np.pi / 2.

        new_lat = np.arcsin(np.sin(lat1) * np.cos(offsetkm/R) + np.cos(lat1) * np.sin(offsetkm/R) * np.cos(bearing))
        new_lon = lon1 + np.arctan2(np.sin(bearing)*np.sin(offsetkm/R)*np.cos(lat1),np.cos(offsetkm/R)-np.sin(lat1)*np.sin(new_lat))
        end_lon = 180.0 * new_lon / np.pi
        end_lat = 180.0 * new_lat / np.pi
        return end_lat, end_lon

   


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
    if Type == 'ALL': P_list = ['FFld','LIDR','Prb1','Prb2','WinS','CoMeT1','CoMeT2','CoMeT3','UAS', 'Ka1','Ka2','NOXP','WSR88D']
    elif Type == 'RADAR': P_list = ['Ka1','Ka2','NOXP','WSR88D']
    elif Type == 'TInsitu': P_list = ['FFld','LIDR','Prb1','Prb2','WinS','CoMeT1','CoMeT2','CoMeT3','UAS']
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
    elif pname == "WTx Mesonet": marker_style, marker_color, line_color, legend_str= r'$\AA$' ,'black','black','WTM Site'   
        #U+1278, #u20a9
    elif pname == 'Ka2': marker_style, marker_color, line_color, legend_str= '8','xkcd:crimson','xkcd:crimson','Ka2'
    elif pname == 'Ka1': marker_style, marker_color, line_color, legend_str= '8','mediumseagreen','mediumseagreen','Ka1'
    elif pname == 'WSR88D': marker_style, marker_color, line_color, legend_str= r'$\Omega$', 'white','black', "WSR88D"

    ##create the legend elements that will be appended to the legend array    
    if pname in pform_names('KA'): #the legend entries for the KA radars
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor='black',markeredgewidth=3,label=legend_str,markerfacecolor=marker_color, markersize=26)
    elif pname == 'WSR88D':
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor=marker_color,markeredgewidth=3,label=legend_str, markersize=26,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')])
    elif pname == 'WTxM':
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor=marker_color,label=legend_str, markersize=26)
    else:
        legend_entry=Line2D([], [], marker=marker_style, markeredgecolor=marker_color,markeredgewidth=3,label=legend_str, markersize=26,path_effects=[PathEffects.withStroke(linewidth=12,foreground='k')])
        
    if print_long== True: print('Made it through pform_attr')
    return marker_style,marker_color,line_color,legend_str,legend_entry

#**************
def read_TInsitu(pname,day,print_long,e_test,tstart=None,tend=None,d_testing=False):
    ''' Reads data files provided on TORUS19 EOL site (Important that datafiles and filenames follow TORUS19 readme)
        ---
        INPUT: filename string (following readme conventions for each platform) 
               tstart,tend: datetime objects
        OUTPUT: dataframe containing all data collected during desired times
    '''
    if pname in pform_names('UNL'): 
        mmfile= glob.glob(temploc+'TORUS_Data/'+day+'/mesonets/UNL/UNL.'+pname+'.*')
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
        mmfile=glob.glob(filesys+'TORUS_Data/'+day+'/mesonets/NSSL/'+pname+'_'+day[2:]+'_QC_met.dat')
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
        return data_nssl, 'NSSL'
    
    # * * * 
    elif pname == 'UAS': 
        if print_long == True: print("no code for UAS yet")
        if d_testing==True: return False
        # + + + + + + + + + + + + ++ + +
	return 'UAS' 

#**************
def radar_dep(pname,scantime,day,print_long,e_test,d_testing=False): 
    ''' Determine if a given radar is deployed and if so assign the correct location values to it.
    '''
    if pname in pform_names('KA'): 
        ## Read in files 
        if pname == 'Ka1': r_testing='ka1' # r_testing is the name of the radar you are testing to see if deployed
        elif pname == 'Ka2': r_testing='ka2'
        #  read in the csv file; if the Ka radar didn't deploy there will be no csv file and the defn will fail
        kadep=pd.read_csv(filesys+'TORUS_Data/'+day+'/radar/TTUKa/csv/'+day+'_deployments_'+r_testing+'.csv')
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
                if p.type == 'NSSL': p.min_max(p_var,mask=True)
                else: p.min_max(p_var)
    
    # * * * 
    elif DType == 'RADAR':
        #this can change for each image (aka from diff scans of the plotting radar) so will be updated for every plot  
        #  ie. are the radars deployed at a given time which (if any) WSR88D fall within plotting domain etc 
        for pname in pform_names('RADAR'):
            if pname in pform_names('KA'):
                data_avail=Platform.test_data(pname,Data)

                if data_avail==True:
                    #  dont repeatedly append to the list if the platform is already included
                    if pname in subset_pnames: pass 
                    else: subset_pnames.append(pname) #append the pname to the subset_pnames list 
                    #  load loc data (and in this case rhi angles if applicable)
                    Data.update({pname: Radar(pname,Data)}) 
                    if print_long==True: print("Data can be read in for platform %s" %(pname))
                else: 
                    if print_long==True: print("No data available to be read in for platform %s" %(pname))

            elif pname in ['NOXP','WSR88D']: print('code not written yet')
    
    # * * * 
    elif DType == 'Plotting_Radar':
        #this oject contains the data relevant to the particular radar data being plotted for each image 
        #  (aka platform name, scantime, the data file etc) this will be updated for each image
        Data.update({'P_Radar':Radar(rname, Rfile=rfile,Swp=swp,Plotting_Radar=True)})

    # * * * 
    elif DType == 'PVar': 
        #Create a Pvar object add it to the data dict and find the global max and min
        #  This object will only be once updated once (aka same for all images for a given run)
        Data.update({p_var: Pvar(p_var)})
        Data[p_var].find_global_max_min(Data)
        Data[p_var].make_dummy_plot()

    # * * * 
    ##Uncomment to check yourself
    # ***************************
    #  print(subset_pnames)
    #  print(Data)
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
                data_avail=read_TInsitu(pname, self.Day, self.Print_long, self.E_test, self.Tstart, self.Tend, d_testing=True)
            elif pname in pform_names('RADAR'):
                data_avail=radar_dep(pname, Data['P_Radar'].time, self.Day, self.Print_long, self.E_test, d_testing=True)
            return data_avail
        except: 
            error_printing(e_test)
            return False

####
class Torus_Insitu(Platform):
    def __init__(self, Name):
        #read in and intilize the dataset; type is the platform type (aka NSSL, UNL, UAS) 
        self.df, self.type= read_TInsitu(Name, self.Day, self.Print_long, self.E_test, self.Tstart, self.Tend)
        self.df.name='{}_df'.format(Name) #assign a name to the pandas dataframe itself
        Platform.__init__(self, Name)

    def min_max(self,p_var,mask=False):
        #Mask the dataset(if applicable) and determine the max and min values of pvar: 
        #  Masking is based of QC flags etc and can be diff for each platform
        if self.type == "NSSL":
            if mask== True: #mask is a bool, sent to False to by default
                self.mask_df= np.ma.masked_where(self.df['qc_flag'].values>0, self.df[p_var].values) #add masked dataset to the object
                self.Min, self.Max = self.mask_df.min(), self.mask_df.max()#use the masked dataset to find max and min of p_var
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
                    
####
class Radar(Platform):
    def __init__(self, Name, Data=None, Rfile= None, Swp=None, Plotting_Radar= False):
        if Plotting_Radar == False:
            self.R_Loc, self.type = radar_dep(Name, Data['P_Radar'].time, self.Day, self.Print_long, self.E_test)
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

#########################################
### set up Pvar class (this is not a subclass of Platform)
class Pvar:
    def __init__(self,p_var):
        self.name=p_var 

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

