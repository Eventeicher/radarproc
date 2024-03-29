#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import needed modules
######################
import matplotlib
#  matplotlib.use('agg')
#  matplotlib.use('GTK3Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import PercentFormatter
import pandas as pd
import numpy as np
import datetime as dt
import csv 
import math 
import pprint
import glob
from datetime import datetime, timedelta
from scipy import stats
from scipy import ndimage, interpolate
from shapely.geometry import LineString, Point, Polygon, MultiPoint
import os
import os.path
import copy
import numpy.ma as ma
import matplotlib.colors as colors
#  from matplotlib.colors import DivergingNorm
from pathlib import Path
import cmocean
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import statistics
import statsmodels.api as sm 

#this is the file with the plotting controls to access any of the vars in that file use config.var
import config as plot_config
from converge_plt import make_con_csv, csvstats_scatterplts, Conv_1dplt, Conv_histograms
from read_pforms import Platform, time_in_range
from radar_defns import get_WSR_from_AWS, read_from_radar_file, radar_fields_prep

def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

def clicker_defn(PLT, AXES):
    pts = [np.nan]
    # Now do a zoom
    while True:
        tellme('Select two corners of zoom, click enter when done')
        pts = plt.ginput(2, timeout=-1)
        if len(pts) < 2:
            break
        (x0, y0), (x1, y1) = pts
        xmin, xmax = sorted([x0, x1])
        ymin, ymax = sorted([y0, y1])
        for a in AXES[::2]: 
            plt.sca(a)
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)

    tellme('Now select surge points,\n click enter when done\n')
    pts = plt.ginput(n=40,timeout=0) # look up ginput docs for some good guidance
    #  plt.waitforbuttonpress()
    return pts

def compute_surgestats(self, day, Data, config, Surge_ID, Surge_Holder, Conv, valid_rays):
    Surgedata=Surge_Holder[Surge_ID]
    #  fields = ['Day','RName', 'Scan_time', 'Surge_ID', 'Surge_Rorigin', 'Pnts_Selected',
            #  'Max_Dist_Offset','AllRay_Count', 'Vray_Count','Straight_Length','length','area', 'slope','Intercept','Max_Vel', 'Min_Vel','Vel_Diff',
            #  'Sigma','Mu','Kappa','Gamma', 'Vray_Sigma', 'Vray_Mu', 'Vray_Kappa','Vray_Gamma', 'deg_thres',
            #  'Conv_Max_Vrays', 'Conv_Max_Allrays', 'Conv_Mean_Vrays', 'Conv_Mean_Allrays', 'Conv_PTP_Vrays', 'Conv_PTP_Allrays',
            #  'Vdiff_convcalc_Max_Vrays', 'Vdiff_convcalc_Max_Allrays', 'Vdiff_convcalc_Mean_Vrays', 'Vdiff_convcalc_Mean_Allrays', 'Vdiff_convcalc_PTP_Vrays', 'Vdiff_convcalc_PTP_Allrays']
    Surge_stats_dict={}
    #  for i in fields:
        #  Surge_stats_dict[i]=None

    # * * *
    #current_surge_df, rname, Surges
    current_surge_df= self.valid_surges[self.valid_surges['Surge_ID'] == Surge_ID]
    rname=current_surge_df.loc[0,'RName']
    Surge_stats_dict['Day']=day
    Surge_stats_dict['RName']=Data['P_Radar'].site_name
    Surge_stats_dict['Surge_ID']=Surge_ID
    Surge_stats_dict['Scan_time']=Data['P_Radar'].Scan_time

    Surge_stats_dict['Surge_Rorigin']=str(rname)+' '+str(current_surge_df.loc[0,'Tilt'])
    Surge_stats_dict['Pnts_Selected']=len(Surgedata['point_labels'])
    Surge_stats_dict['Max_Dist_Offset']=np.max(config.Surge_controls['Surges']['offset_dist'])
    Surge_stats_dict['AllRay_Count']=np.max(Surgedata['Center_pnts']['Index']['Ray'])-np.min(Surgedata['Center_pnts']['Index']['Ray'])
    Surge_stats_dict['VRay_Count']=(len(Surgedata['Surge_Subset']['Valid_Rays']))#, label=self.Surges[Surge_ID]['point_labels'][pnt])

    init_point = Point(Surgedata['Center_pnts']['x'][0], Surgedata['Center_pnts']['y'][0])
    final_point = Point(Surgedata['Center_pnts']['x'][-1], Surgedata['Center_pnts']['y'][-1])
    straight_dist=init_point.distance(final_point)
    Surge_stats_dict['Straight_Length']=round(straight_dist)
    Surge_stats_dict['length']=round(Surgedata['Center_pnts']['line_object'].length)
    Surge_stats_dict['area']=round(Surgedata['Polygon'].area)

    if Surgedata['line']['slope'] != None:  
        Surge_stats_dict['slope']=round(Surgedata['line']['slope'],3)
        Surge_stats_dict['Intercept']=round(Surgedata['line']['intercept'])
        Surge_stats_dict['Max_Vel']=round(Surgedata['Surge_Subset']['Max']['value'])
        Surge_stats_dict['Min_Vel']=round(Surgedata['Surge_Subset']['Min']['value'])
        Surge_stats_dict['Vel_Diff']=abs(Surge_stats_dict['Max_Vel']-Surge_stats_dict['Min_Vel'])

    flat_Conv=np.ndarray.flatten(Conv)
    Conv_nonans=flat_Conv[~np.isnan(flat_Conv)]

    Surge_stats_dict['Sigma']=np.std(Conv_nonans)
    Surge_stats_dict['Mu']=np.mean(Conv_nonans)
    Surge_stats_dict['Kappa']=stats.kurtosis(Conv_nonans)
    Surge_stats_dict['Gamma']=stats.skew(Conv_nonans)
    
    if config.Conv_controls['Valid_Rays'] == True:
        Surge_stats_dict['deg_thres']=str(config.Surge_controls['Surges']['ray_selection_method'])+' '+str(config.Surge_controls['Surges']['ray_selection_thresh'])
        Vr_Conv=Conv
        for ray in range(np.shape(Vr_Conv)[0]):
            if valid_rays[ray]==False:
                Vr_Conv[ray,:]=np.nan
        flat_Vr_Conv=np.ndarray.flatten(Vr_Conv)
        Vr_Conv_nonans=flat_Vr_Conv[~np.isnan(flat_Vr_Conv)]
        
        if len(Vr_Conv_nonans)!=0:
            Surge_stats_dict['Vray_Sigma']=np.std(Vr_Conv_nonans)
            Surge_stats_dict['Vray_Mu']=np.mean(Vr_Conv_nonans)
            Surge_stats_dict['Vray_Kappa']=stats.kurtosis(Vr_Conv_nonans)
            Surge_stats_dict['Vray_Gamma']=stats.skew(Vr_Conv_nonans)

            Surge_stats_dict['Sigma_Diff']=abs(Surge_stats_dict['Sigma']-Surge_stats_dict['Vray_Sigma'])
            Surge_stats_dict['Mu_Diff']=abs(Surge_stats_dict['Mu']-Surge_stats_dict['Vray_Mu'])
            Surge_stats_dict['Kappa_Diff']=abs(Surge_stats_dict['Kappa']-Surge_stats_dict['Vray_Kappa'])
            Surge_stats_dict['Gamma_Diff']=abs(Surge_stats_dict['Gamma']-Surge_stats_dict['Vray_Gamma'])

    Surge_stats_dict['Conv_Max_Vrays']=Surgedata['Diff_Calcs']['Conv']['Max']['Vrays']
    Surge_stats_dict['Conv_Max_Allrays']=Surgedata['Diff_Calcs']['Conv']['Max']['Allrays']
    Surge_stats_dict['Conv_Mean_Vrays']=Surgedata['Diff_Calcs']['Conv']['Mean']['Vrays']
    Surge_stats_dict['Conv_Mean_Allrays']=Surgedata['Diff_Calcs']['Conv']['Mean']['Allrays']

    Surge_stats_dict['Vdiff_convcalc_Max_Vrays']=Surgedata['Diff_Calcs']['Vdiff_convcalc']['Max']['Vrays']
    Surge_stats_dict['Vdiff_convcalc_Max_Allrays']=Surgedata['Diff_Calcs']['Vdiff_convcalc']['Max']['Allrays']
    Surge_stats_dict['Vdiff_convcalc_Mean_Vrays']=Surgedata['Diff_Calcs']['Vdiff_convcalc']['Mean']['Vrays']
    Surge_stats_dict['Vdiff_convcalc_Mean_Allrays']=Surgedata['Diff_Calcs']['Vdiff_convcalc']['Mean']['Allrays']
    # * * * *
    if len(Surge_stats_dict['Conv_Max_Vrays']) != 0:
        #  print(Surge_stats_dict['Conv_Max_Vrays'])
        #  print(type(Surge_stats_dict['Conv_Max_Vrays']))
        #  print('7777777777777777777')
        #  print(np.mean(Surge_stats_dict['Conv_Max_Vrays']))
        #  print(type(np.mean(Surge_stats_dict['Conv_Max_Vrays'])))
        #  print('8888888887777777777')
        Surge_stats_dict['Conv_Max_Vrays_Mean']=np.mean(Surge_stats_dict['Conv_Max_Vrays'])
        Surge_stats_dict['Conv_Mean_Vrays_Mean']=np.mean(Surge_stats_dict['Conv_Mean_Vrays'])
        Surge_stats_dict['Conv_Max_Vrays_Median']=np.median(Surge_stats_dict['Conv_Max_Vrays'])
        Surge_stats_dict['Conv_Mean_Vrays_Median']=np.median(Surge_stats_dict['Conv_Mean_Vrays'])
        Surge_stats_dict['Conv_Max_Vrays_Min']=np.min(Surge_stats_dict['Conv_Max_Vrays'])
        Surge_stats_dict['Conv_Mean_Vrays_Min']=np.min(Surge_stats_dict['Conv_Mean_Vrays'])
        Surge_stats_dict['Conv_Max_Vrays_Max']=np.max(Surge_stats_dict['Conv_Max_Vrays'])
        Surge_stats_dict['Conv_Mean_Vrays_Max']=np.max(Surge_stats_dict['Conv_Mean_Vrays'])
        Surge_stats_dict['Vdiff_convcalc_Max_Vrays_Mean']=np.mean(Surge_stats_dict['Vdiff_convcalc_Max_Vrays'])
        Surge_stats_dict['Vdiff_convcalc_Mean_Vrays_Mean']=np.mean(Surge_stats_dict['Vdiff_convcalc_Mean_Vrays'])
        Surge_stats_dict['Vdiff_convcalc_Max_Vrays_Median']=np.median(Surge_stats_dict['Vdiff_convcalc_Max_Vrays'])
        Surge_stats_dict['Vdiff_convcalc_Mean_Vrays_Median']=np.median(Surge_stats_dict['Vdiff_convcalc_Mean_Vrays'])
        Surge_stats_dict['Vdiff_convcalc_Max_Vrays_Max']=np.max(Surge_stats_dict['Vdiff_convcalc_Max_Vrays'])
        Surge_stats_dict['Vdiff_convcalc_Mean_Vrays_Max']=np.max(Surge_stats_dict['Vdiff_convcalc_Mean_Vrays'])
        Surge_stats_dict['Vdiff_convcalc_Max_Vrays_Min']=np.min(Surge_stats_dict['Vdiff_convcalc_Max_Vrays'])
        Surge_stats_dict['Vdiff_convcalc_Mean_Vrays_Min']=np.min(Surge_stats_dict['Vdiff_convcalc_Mean_Vrays'])

    Surge_stats_dict['Conv_Max_Allrays_Mean']=np.mean(Surge_stats_dict['Conv_Max_Allrays'])
    Surge_stats_dict['Conv_Mean_Allrays_Mean']=np.mean(Surge_stats_dict['Conv_Mean_Allrays'])

    Surge_stats_dict['Conv_Max_Allrays_Median']=np.median(Surge_stats_dict['Conv_Max_Allrays'])
    Surge_stats_dict['Conv_Mean_Allrays_Median']=np.median(Surge_stats_dict['Conv_Mean_Allrays'])

    Surge_stats_dict['Conv_Max_Allrays_Max']=np.max(Surge_stats_dict['Conv_Max_Allrays'])
    Surge_stats_dict['Conv_Mean_Allrays_Max']=np.max(Surge_stats_dict['Conv_Mean_Allrays'])

    Surge_stats_dict['Conv_Max_Allrays_Min']=np.min(Surge_stats_dict['Conv_Max_Allrays'])
    Surge_stats_dict['Conv_Mean_Allrays_Min']=np.min(Surge_stats_dict['Conv_Mean_Allrays'])
    
    Surge_stats_dict['Vdiff_convcalc_Max_Allrays_Mean']=np.mean(Surge_stats_dict['Vdiff_convcalc_Max_Allrays'])
    Surge_stats_dict['Vdiff_convcalc_Mean_Allrays_Mean']=np.mean(Surge_stats_dict['Vdiff_convcalc_Mean_Allrays'])

    Surge_stats_dict['Vdiff_convcalc_Max_Allrays_Median']=np.median(Surge_stats_dict['Vdiff_convcalc_Max_Allrays'])
    Surge_stats_dict['Vdiff_convcalc_Mean_Allrays_Median']=np.median(Surge_stats_dict['Vdiff_convcalc_Mean_Allrays'])

    Surge_stats_dict['Vdiff_convcalc_Max_Allrays_Max']=np.max(Surge_stats_dict['Vdiff_convcalc_Max_Allrays'])
    Surge_stats_dict['Vdiff_convcalc_Mean_Allrays_Max']=np.max(Surge_stats_dict['Vdiff_convcalc_Mean_Allrays'])

    Surge_stats_dict['Vdiff_convcalc_Max_Allrays_Min']=np.min(Surge_stats_dict['Vdiff_convcalc_Max_Allrays'])
    Surge_stats_dict['Vdiff_convcalc_Mean_Allrays_Min']=np.min(Surge_stats_dict['Vdiff_convcalc_Mean_Allrays'])


    ####### 
    print('HHLLLLLLLLL')
    print(Surge_stats_dict.keys())
    print('JJJJKKJJJJ')
    make_surgestats_csv(config, Data, Surge_stats_dict.keys(), Surge_stats_dict)
    pprint.pprint(Surge_stats_dict)
    return Surge_stats_dict

def make_surgestats_csv(config, Data, fields, Surge_stats_dict):
    csv_name= config.g_TORUS_directory+'Surge_stats.csv'
    Does_csv_exist=os.path.isfile(csv_name)
    if Does_csv_exist == False: 
        # making csv file if it does not exist 
        with open(csv_name, 'w') as csvfile: 
            # field names 
            
            # creating a csv writer object 
            csvwriter = csv.writer(csvfile) 
            # writing the fields 
            csvwriter.writerow(fields) 

    ###
    # Now append the new row of data points
    with open(csv_name, 'a') as csvfile: 
        writer= csv.DictWriter(csvfile, fieldnames=fields)
        #  for data in Surge_stats_dict:
        writer.writerow(Surge_stats_dict)


def make_surgemeso_obj_csv(config, Data, pts, surge_name, tilt, filename, Type, Does_csv_exist):
    if Type == 'Surge':
        ID_name = 'Surge_ID'
    elif Type == 'Meso':
        ID_name = 'Meso_ID'

    sweep = Data['P_Radar'].swp
    if config.Radar_Plot_Type == 'WSR_Plotting':
        sweep= sweep[1]
    gate_lat, gate_lon, _ = Data['P_Radar'].rfile.get_gate_lat_lon_alt(sweep)
    gate_x, gate_y, _ = Data['P_Radar'].rfile.get_gate_x_y_z(sweep)

    if Does_csv_exist == False: 
        # making csv file if it does not exist 
        with open(filename, 'w') as csvfile: 
            # field names 
            fields = ['RName', 'Tilt', 'Scan_time', ID_name, 'point_x', 'point_y', 'point_lat', 'point_lon'] 
            
            # creating a csv writer object 
            csvwriter = csv.writer(csvfile) 
            # writing the fields 
            csvwriter.writerow(fields) 
    ###
    # Now append the new row of data points
    with open(filename, 'a') as csvfile: 
        # writing the data rows 
        csvwriter = csv.writer(csvfile) 
        for i in np.arange(0, len(pts)):
            point=pts[i]
        
            if Data['P_Radar'].name == 'WSR88D':
                Rname= Data['P_Radar'].site_name
            else: 
                Rname= Data['P_Radar'].name
            
            cray_ind, crange_ind, _ = return_pnt_index(Data, gate_x, gate_y, sweep, given_x_lon=point[0], given_y_lat=point[1])
            p_lat, p_lon = gate_lat[cray_ind, crange_ind], gate_lon[cray_ind, crange_ind]
            rows= [[Rname, tilt, Data['P_Radar'].Scan_time, surge_name, point[0],point[1], p_lat, p_lon]]
            csvwriter.writerows(rows)

def return_pnt_index(Data, gate_x_lon, gate_y_lat, sweep, Buffer_dist=None, surge_pd=None, given_x_lon=None, given_y_lat=None, point=None, XY_or_LatLon='xy'):
    def find_index_of_nearest_xy(y_array, x_array, y_point, x_point):
        distance = (y_array-y_point)**2 + (x_array-x_point)**2
        idy,idx = np.where(distance==distance.min())
        return idy[0],idx[0]
    # * * * *
    if isinstance(surge_pd, pd.DataFrame):
        if XY_or_LatLon == 'xy':
            pnt_X_Lon, pnt_Y_Lat = surge_pd.point_x[point], surge_pd.point_y[point]
        if XY_or_LatLon == 'latlon':
            pnt_X_Lon, pnt_Y_Lat = surge_pd.point_lon[point], surge_pd.point_lat[point]
    else:
        pnt_X_Lon, pnt_Y_Lat = given_x_lon, given_y_lat 
    ray_ind, range_ind = find_index_of_nearest_xy(gate_y_lat, gate_x_lon, pnt_Y_Lat, pnt_X_Lon)

    if Buffer_dist != None:
        bin_dist= Data['P_Radar'].rfile.range['meters_between_gates']
        extra_bins = round(Buffer_dist/bin_dist)
    else: extra_bins = None
    #  self.display.plot_line_xy(surge_pd.point_x, surge_pd.point_y)
    #  pnt=points.split("), (")
    return ray_ind, range_ind, extra_bins

def csv_reader(config, Data, day, tilt, Type):
    def subset_df(df, Data, ID, tilt, Type):
        Crossing_Radars = config.Surge_controls['Feature_IDing']['Existing_pnts']['cross_rpforms']['allow']
        def no_valid_object_checker(df):
            if len(df) == 0: valid_object= False
            else:  valid_object= True
            return valid_object
        def nearest(items, pivot):
            return min(items, key=lambda x: abs(x - pivot))
        # # # # 
        if Type == 'Surge':  object_ID = Type+'_ID'
        elif Type == 'Meso':  object_ID = Type+'_ID'
        ID_sub = df[df[object_ID] == ID]
        remaining_objects= no_valid_object_checker(ID_sub)

        if remaining_objects == True:
            if Crossing_Radars == False: valid_time = str(Data['P_Radar'].Scan_time)
            else:
                try:
                    times_list = [dt.datetime.strptime(time, '%Y-%m-%d %H:%M:%S') for time in df.Scan_time.unique()]
                except:
                    times_list = [dt.datetime.strptime(time, '%Y-%m-%d %H:%M:%S.%f') for time in df.Scan_time.unique()]
                closest_time=nearest(times_list, Data['P_Radar'].Scan_time)
                max_time_range = Data['P_Radar'].Scan_time + dt.timedelta(minutes=config.Surge_controls['Feature_IDing']['Existing_pnts']['cross_rpforms']['within_time'])
                min_time_range = Data['P_Radar'].Scan_time - dt.timedelta(minutes=config.Surge_controls['Feature_IDing']['Existing_pnts']['cross_rpforms']['within_time'])
                if time_in_range(min_time_range, max_time_range, closest_time): valid_time = str(closest_time) 
                else: valid_time = None

            time_sub=ID_sub[ID_sub['Scan_time'] == valid_time]
            remaining_objects= no_valid_object_checker(time_sub)
        
        if remaining_objects == True:
            if Crossing_Radars == False: valid_radar = [Data['P_Radar'].name]
            else: 
                valid_radar = time_sub.RName.unique()

            radar_sub=time_sub[time_sub['RName'].isin(valid_radar)]
            remaining_objects= no_valid_object_checker(radar_sub)

        if remaining_objects == True:
            if Crossing_Radars == False: valid_tilt = tilt 
            else: 
                tilt_allowance = config.Surge_controls['Feature_IDing']['Existing_pnts']['cross_rpforms']['tilt_allowance']
                if tilt_allowance == None: valid_tilt= tilt 
                else: 
                    possible_tilts = radar_sub.Tilt.unique()
                    closest_tilt = nearest(possible_tilts, tilt)
                    max_tilt, min_tilt = (tilt + tilt_allowance), (tilt - tilt_allowance)
                    if (closest_tilt <= max_tilt) and (closest_tilt >= min_tilt): valid_tilt = closest_tilt
                    else:  valid_tilt = None

            tilt_sub=radar_sub[radar_sub['Tilt'] == valid_tilt]
            remaining_objects= no_valid_object_checker(tilt_sub)
        
        if remaining_objects == True: final_sub = tilt_sub.reset_index()
        else: final_sub = None
        return final_sub
        
    # ****
    try: 
        if Type == 'Surge': object_type = 'surge'
        elif Type == 'Meso': object_type = 'meso'
        df=pd.read_csv(config.g_TORUS_directory+day+'/data/'+day+'_'+object_type+'_pnts.csv')
    except: 
        print('did not find file: '+ config.g_TORUS_directory+day+'/data/'+day+'_'+object_type+'_pnts.csv' )
        print('NO '+ Type+ ' READ IN FOR THIS DAY')
        num_of_objects= 0
        return num_of_objects, None , None

    if Type == 'Surge': objects=df.Surge_ID.unique()
    elif Type == 'Meso': objects=df.Meso_ID.unique()

    num_of_objects, object_ids, objects_attime_df = 0, [], pd.DataFrame() 
    for j in objects: 
        single_object_df=subset_df(df, Data, j, tilt, Type)
        if isinstance(single_object_df, pd.DataFrame):
            objects_attime_df=objects_attime_df.append(single_object_df)
            num_of_objects=num_of_objects+1
            object_ids.append(j)
    return num_of_objects, object_ids, objects_attime_df  



def meso_radarbins(self, day, config, Data, sweep):
    if config.Radar_Plot_Type=='WSR_Plotting': #comeback 
        sweep=sweep[1]
    Meso_Holder={}

    ######################################3
    plt_gate_x, plt_gate_y, _ = Data['P_Radar'].rfile.get_gate_x_y_z(sweep)
    plt_gate_lat, plt_gate_lon, _ = Data['P_Radar'].rfile.get_gate_lat_lon_alt(sweep)

    if Data['P_Radar'].site_name in self.valid_mesos.RName.unique()[0]:
        origR_gate_x, origR_gate_y = plt_gate_x, plt_gate_y
        origR_gate_lat, origR_gate_lon = plt_gate_lat, plt_gate_lon
        r_mom_data=Data['P_Radar'].rfile.fields['Meso_azi']['data']
        sweep_startidx = np.int64(Data['P_Radar'].rfile.sweep_start_ray_index['data'][sweep])
        sweep_endidx = np.int64(Data['P_Radar'].rfile.sweep_end_ray_index['data'][sweep])
    else:
        origR_RNAME=self.valid_mesos.RName.unique()[0]
        time=self.valid_mesos.Scan_time.unique()[0]
        print(time)
        TIME=dt.datetime.strptime(time, '%Y-%m-%d %H:%M:%S.%f') 
        if config.radar_controls['Use_downloaded_files']== True:
            path = config.g_TORUS_directory + day+'/data/radar/Nexrad/Nexrad_files/dealiased_'+origR_RNAME+'*'+time[11:13]+time[14:16]+'*'
            print(path)
            r_files_path = sorted(glob.glob(path))
            print(r_files_path)
        else:
            radar_file = get_WSR_from_AWS(config, day, (TIME-timedelta(minutes=1)), (TIME+timedelta(minutes=1)), 'KLBB')
        radar = read_from_radar_file(config, radar_file[0], True)
        radar = radar_fields_prep(config, radar, 'WSR', 1, moment='Meso_azi')
        Meso_Holder.update({'RADAR':{'name':origR_RNAME, 'Scan_time':TIME, 'rfile':radar, 'swp':[0,1]}})

        origR_gate_x, origR_gate_y, _ = radar.get_gate_x_y_z(1)
        origR_gate_lat, origR_gate_lon, _ = radar.get_gate_lat_lon_alt(1)
        r_mom_data=radar.fields['Meso_azi']['data']
        sweep_startidx = np.int64(radar.sweep_start_ray_index['data'][1])
        sweep_endidx = np.int64(radar.sweep_end_ray_index['data'][1])
    ######################################3
    for Meso_ID in self.meso_ids[:]:
        zoom_meso_df= self.valid_mesos[self.valid_mesos['Meso_ID'] == Meso_ID]
        if Data['P_Radar'].site_name in zoom_meso_df.RName.unique():
            pass
        else: 
            for i in zoom_meso_df.index.tolist():
                ray_ind, range_ind, extra_bins = return_pnt_index(Data, plt_gate_lon, plt_gate_lat, sweep, surge_pd=zoom_meso_df, 
                                                                  given_x_lon=zoom_meso_df.loc[i,'point_lon'],
                                                                  given_y_lat=zoom_meso_df.loc[i,'point_lat'], XY_or_LatLon='latlon', point=i)
                zoom_meso_df.loc[i, 'point_x']= plt_gate_x[ray_ind, range_ind]
                zoom_meso_df.loc[i, 'point_y']= plt_gate_y[ray_ind, range_ind]
 
        # * * *
        pnt_name = []
        rays_ind, ranges_ind = [], []
        r_bin_sub, r_data_sub = [], []
        for i in zoom_meso_df.index.tolist():
            pnt_name.append(i)
            Meso_Holder.update({Meso_ID:{'point_labels': pnt_name}})  

            #############
            #determine the index positions for each point along the surge (and for each offset point)
            ray_ind, range_ind, extra_bins = return_pnt_index(Data, plt_gate_x, plt_gate_y, sweep, surge_pd=zoom_meso_df, point=i)
            rays_ind.append(ray_ind)
            ranges_ind.append(range_ind)
            Ray_Index, Range_Index= ray_ind, range_ind
            plt_xpos, plt_ypos = plt_gate_x[Ray_Index , Range_Index], plt_gate_y[Ray_Index , Range_Index]
                   
            point_entry = {'x': plt_xpos, 
                           'y': plt_ypos, 
                           'Point_object': Point(plt_xpos, plt_ypos),
                           'Ring':Point(plt_xpos, plt_ypos).buffer(config.Surge_controls['Mesos']['radius'])
                           }
            Meso_Holder.update({Meso_ID:{'Center_pnt': point_entry}})  
            
            #############
            square_bounds=point_entry['Ring'].bounds
            #  mask=np.ma.getmask(origR_gate_x)
            mask=np.ma.masked_array(origR_gate_x)
            for ray in range(origR_gate_x.shape[0]): #iterate over each ray
                for rbin in range(origR_gate_x.shape[1]): #iterate over each bin 
                    mask_entered= False 
                    if (origR_gate_x[ray,rbin] >= square_bounds[0]) and (origR_gate_x[ray,rbin] <= square_bounds[2]):
                        if (origR_gate_y[ray,rbin] >= square_bounds[1]) and (origR_gate_y[ray,rbin] <= square_bounds[3]):
                            #  (x - center_x)^2 + (y - center_y)^2 < radius^2
                            A=(origR_gate_x[ray,rbin] - plt_xpos)**2 + (origR_gate_y[ray,rbin] - plt_ypos)**2 
                            if A <= config.Surge_controls['Mesos']['radius']**2:
                               mask[ray,rbin] = False
                               mask_entered= True

                    if mask_entered == False:
                        mask[ray,rbin] =True 

            masked_subset_rmom= np.where(mask== False, r_mom_data[sweep_startidx:sweep_endidx+1], np.nan)
            final_subset_rmom= np.where(masked_subset_rmom== -9999.0, np.nan, masked_subset_rmom)
            
            maxindexlist=np.where(final_subset_rmom == np.nanmax(final_subset_rmom))
            minindexlist=np.where(final_subset_rmom == np.nanmin(final_subset_rmom))
            xmax, ymax, xmin, ymin = [],[],[],[]
            if Data['P_Radar'].site_name in zoom_meso_df.RName.unique():
                for pnt in range(len(maxindexlist[0])):
                    xmax.append(origR_gate_x[maxindexlist[0][pnt], maxindexlist[1][pnt]])
                    ymax.append(origR_gate_y[maxindexlist[0][pnt], maxindexlist[1][pnt]])
                for pnt in range(len(minindexlist[0])):
                    xmin.append(origR_gate_x[minindexlist[0][pnt], minindexlist[1][pnt]])
                    ymin.append(origR_gate_y[minindexlist[0][pnt], minindexlist[1][pnt]])
            else: 
                for pnt in range(len(maxindexlist[0])):
                    origR_lat=origR_gate_lat[maxindexlist[0][pnt], maxindexlist[1][pnt]]
                    origR_lon=origR_gate_lon[maxindexlist[0][pnt], maxindexlist[1][pnt]]
                    ray_ind, range_ind, extra_bins = return_pnt_index(Data, plt_gate_lon, plt_gate_lat, sweep, surge_pd=zoom_meso_df, 
                                                                  given_x_lon=origR_lon, given_y_lat=origR_lat, XY_or_LatLon='latlon', point=i)
                    xmax.append(plt_gate_x[ray_ind, range_ind])
                    ymax.append(plt_gate_y[ray_ind, range_ind])
                
                for pnt in range(len(minindexlist[0])):
                    origR_lat=origR_gate_lat[minindexlist[0][pnt], minindexlist[1][pnt]]
                    origR_lon=origR_gate_lon[minindexlist[0][pnt], minindexlist[1][pnt]]
                    ray_ind, range_ind, extra_bins = return_pnt_index(Data, plt_gate_lon, plt_gate_lat, sweep, surge_pd=zoom_meso_df, 
                                                                  given_x_lon=origR_lon, given_y_lat=origR_lat, XY_or_LatLon='latlon', point=i)
                    xmin.append(plt_gate_x[ray_ind, range_ind])
                    ymin.append(plt_gate_y[ray_ind, range_ind])

            Meso_Holder.update({'Max':{'x':xmax,'y':ymax, 'value':np.nanmax(final_subset_rmom)}})
            Meso_Holder.update({'Min':{'x':xmin,'y':ymin, 'value':np.nanmin(final_subset_rmom)}})

            #  Meso_subset_rmom= np.zeros(r_mom_ofintrest_data.shape)
            #  Meso_subset_rmom[sweep_startidx:sweep_endidx+1] = masked_subset_orig_rmom
            #  Data['P_Radar'].rfile.add_field('Meso_test', {'data': Meso_subset_rmom}, replace_existing=True)
            #############

    #  pprint.pprint(Meso_Holder)
    return Meso_Holder

def surge_radarbins(self, config, Data, sweep, day):
    def angle3pt(a, b, c):
        """Counterclockwise angle in degrees by turning from a to c around b
            Returns a float between 0.0 and 360.0"""
        #  if isinstance(Plus_Intersections, Point):
        if isinstance(c, Point):
            C_x, C_y = c.x, c.y
        elif isinstance(c, MultiPoint):
            c= c[0]
            C_x, C_y = c.x, c.y
        else: 
            C_x, C_y = c[0], c[1]
        ang = math.degrees(math.atan2(C_y-b[1], C_x-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
        if b in (a, c):
            raise ValueError("Undefined angle, two identical points", (a, b, c))
        return ang + 360 if ang < 0 else ang
    def pnt_intersection(b1, b2, m1, m2):
        xi = (b1-b2) / (m2-m1)
        yi = m1*xi + b1
        return xi, yi
    

    if config.Radar_Plot_Type == 'WSR_Plotting':
        sweep=sweep[1]
    gate_x, gate_y, _ = Data['P_Radar'].rfile.get_gate_x_y_z(sweep)
    gate_lat, gate_lon, _ = Data['P_Radar'].rfile.get_gate_lat_lon_alt(sweep)
    range_bins_alongray = Data['P_Radar'].rfile.range['data']
    r_mom_ofintrest_data=Data['P_Radar'].rfile.fields[config.lineplt_control['H']['var']]['data']
    sweep_startidx = np.int64(Data['P_Radar'].rfile.sweep_start_ray_index['data'][sweep])
    sweep_endidx = np.int64(Data['P_Radar'].rfile.sweep_end_ray_index['data'][sweep])
    subset_orig_rmom=r_mom_ofintrest_data[sweep_startidx:sweep_endidx+1]

    Surge_Holder={}
    save_Surge_subset=0
    for Surge_ID in self.surge_ids[:]:
        zoom_surge_df= self.valid_surges[self.valid_surges['Surge_ID'] == Surge_ID]
        if Data['P_Radar'].name in zoom_surge_df.RName.unique():
            pass
        else: 
            for i in zoom_surge_df.index.tolist():
                ray_ind, range_ind, extra_bins = return_pnt_index(Data, gate_lon, gate_lat, sweep, surge_pd=zoom_surge_df, 
                                                                  given_x_lon=zoom_surge_df.loc[i,'point_lon'],
                                                                  given_y_lat=zoom_surge_df.loc[i,'point_lat'], XY_or_LatLon='latlon', point=i)
                zoom_surge_df.loc[i, 'point_x']= gate_x[ray_ind, range_ind]
                zoom_surge_df.loc[i, 'point_y']= gate_y[ray_ind, range_ind]
 
        # * * *
        pnt_name = []
        rays_ind, ranges_ind = [], []
        r_bin_sub, r_data_sub = [], []
        for i in zoom_surge_df.index.tolist():
            pnt_name.append(i)
            #determine the index positions for each point along the surge (and for each offset point)
            offset_bins=[]
            for offset_dist in config.Surge_controls['Surges']['offset_dist']:
                ray_ind, range_ind, extra_bins = return_pnt_index(Data, gate_x, gate_y, sweep, surge_pd=zoom_surge_df, Buffer_dist=offset_dist, point=i)
                offset_bins.append(extra_bins)
            rays_ind.append(ray_ind)
            ranges_ind.append(range_ind)
        Surge_Holder.update({Surge_ID:{'point_labels': pnt_name}})  
        Surge_Holder[Surge_ID].update({'Offset_pnts': 
                                            {'Max': {'offset': np.max(offset_bins),'dist': np.max(config.Surge_controls['Surges']['offset_dist'])},
                                             'Min': {'offset': np.min(offset_bins),'dist': np.min(config.Surge_controls['Surges']['offset_dist'])}
                                             }
                                        })

           
        #store the relevant data for each surge (and surge offset) point in a dict 
        for o in range(len(config.Surge_controls['Surges']['offset_dist'])):
            Both_Dir={}
            for offset_dir in ['Center','Plus', 'Minus']: 
                if offset_dir == 'Center':  extra_bins= 0
                elif offset_dir == 'Plus':  extra_bins= offset_bins[o]
                elif offset_dir == 'Minus': extra_bins= offset_bins[o] * -1
                
                r_bin, r_data, xpos, ypos = [], [], [], []
                rays, ranges = [], []
                for pnt in range(len(pnt_name)):
                    Ray_Index, Range_Index= rays_ind[pnt], ranges_ind[pnt]+extra_bins  
                    rays.append(Ray_Index)
                    ranges.append(Range_Index)

                    if Range_Index<= len(range_bins_alongray):
                        r_bin.append(range_bins_alongray[Range_Index])
                        r_data.append(r_mom_ofintrest_data[Ray_Index , Range_Index])
                        xpos.append(gate_x[Ray_Index , Range_Index])
                        ypos.append(gate_y[Ray_Index , Range_Index])
                   
                comb_xy=list(zip(xpos, ypos))
                point_entry = {'x': xpos, 
                               'y': ypos, 
                               'range_value': r_bin, 
                               'radar_value': r_data,
                               'line_object': LineString(comb_xy), 
                               'Index':{'Ray': rays, 'Range': ranges}
                               }

                if offset_dir == 'Center': 
                    Surge_Holder[Surge_ID].update({'Center_pnts': point_entry})
                else:
                    Both_Dir.update({offset_dir:point_entry})
                    Surge_Holder[Surge_ID]['Offset_pnts'].update({config.Surge_controls['Surges']['offset_dist'][o]: Both_Dir}) 
                

        #add the subset for the surge area (max offset) 
        Max_obins=np.max(offset_bins)
        for pnt in range(len(pnt_name)):
            #range_bins within surge for a given ray 
            r_bin_sub.append(range_bins_alongray[ranges_ind[pnt]-Max_obins:ranges_ind[pnt]+Max_obins])
            #data within surge singleray
            r_data_sub.append(r_mom_ofintrest_data[rays_ind[pnt],ranges_ind[pnt]-Max_obins:ranges_ind[pnt]+Max_obins])
        Surge_Holder[Surge_ID].update({'Subsets':
                                             {'Radar_data':np.asarray(r_data_sub), 'Range_bins':np.asarray(r_bin_sub)}
                                         }) 

        #########
        #add Polygon for the full (maxoffset) area 
        Max_odist=np.max(config.Surge_controls['Surges']['offset_dist'])
        Xplus=Surge_Holder[Surge_ID]['Offset_pnts'][Max_odist]['Plus']['x']
        Yplus=Surge_Holder[Surge_ID]['Offset_pnts'][Max_odist]['Plus']['y']
        Xminus=Surge_Holder[Surge_ID]['Offset_pnts'][Max_odist]['Minus']['x']
        Yminus=Surge_Holder[Surge_ID]['Offset_pnts'][Max_odist]['Minus']['y']
        
        comb_xy_minus, comb_xy_plus = list(zip(Xminus, Yminus)), list(zip(Xplus, Yplus))
        poly=Polygon(comb_xy_plus+comb_xy_minus)
        Surge_Holder[Surge_ID].update({'Polygon': poly}) 

        # * * *
        if 'zoom' in config.r_mom: 
            c_lon, c_lat = zoom_surge_df.point_lon.mean(), zoom_surge_df.point_lat.mean()

            if self.config.radar_controls['offsetkm']['Zoom'] == 'Fitted':
                surge_extents=Surge_Holder[Surge_ID]['Polygon'].envelope.bounds
                x_dist=(surge_extents[2]-surge_extents[0])/2000
                y_dist=(surge_extents[3]-surge_extents[1])/2000
                if x_dist > y_dist: o_km=x_dist+.4 #.1 is to add an additional buffer
                elif y_dist >= x_dist: o_km=y_dist+.4
            else:
                o_km=self.config.radar_controls['offsetkm']['Zoom']

            ZOOM_Dom = Platform.getLocation('ClickPoint', scan_time=Data['P_Radar'].Scan_time, 
                                              offsetkm=o_km , pnt_x=c_lon, pnt_y=c_lat)
        else: 
            ZOOM_Dom = None

        # * * *
        points = [Surge_Holder[Surge_ID]['Center_pnts']['line_object'].interpolate(dist) for dist in 
                  np.arange(0, Surge_Holder[Surge_ID]['Center_pnts']['line_object'].length, 0.5)]
        if len(points) != 0: 
            x, y = zip(*[(pt.x, pt.y) for pt in points])
            surge_slope, surge_yintercept, r_value, p_value, std_err = stats.linregress(x,y)

            x_vals = (np.min(gate_x), np.max(gate_x))
            x_vals= np.asarray(x_vals)
            y_vals = surge_yintercept + surge_slope * x_vals
            surge_linepoint= (x_vals[-1], y_vals[-1]) #a point along the surge line
        else: 
            surge_slope, surge_yintercept=None, None


        # * * *
        #define the indices for the required sweep
        ray_selection= np.zeros(gate_x.shape)
        rayangle_tosurge=[]
        store_rays_wthin_range=[]
        Surge_Area_Rdata=np.zeros(gate_x.shape)
        Max_dist = Surge_Holder[Surge_ID]['Offset_pnts']['Max']['dist']
        for c_ray in range(ray_selection.shape[0]): #iterate over each ray
            if c_ray in np.arange(np.min(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray']),
                                  np.max(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray'])):
                ray_line_object=LineString(list(zip(gate_x[c_ray], gate_y[c_ray])))

                #  for dist in config.Surge_controls['Surges']['offset_dist']:
                    #  Plus_Intersections=Surge_Holder[Surge_ID]['Offset_pnts'][dist]['Plus']['line_object'].intersection(ray_line_object)
                    #  Minus_Intersections=Surge_Holder[Surge_ID]['Offset_pnts'][dist]['Minus']['line_object'].intersection(ray_line_object)

                Plus_Intersections=Surge_Holder[Surge_ID]['Offset_pnts'][Max_dist]['Plus']['line_object'].intersection(ray_line_object)
                Minus_Intersections=Surge_Holder[Surge_ID]['Offset_pnts'][Max_dist]['Minus']['line_object'].intersection(ray_line_object)

                if isinstance(Plus_Intersections, MultiPoint): Plus_Intx=Plus_Intersections[0]
                elif isinstance(Plus_Intersections, Point): Plus_Intx=Plus_Intersections
                if isinstance(Minus_Intersections, MultiPoint): Minus_Intx=Minus_Intersections[0]
                elif isinstance(Minus_Intersections, Point): Minus_Intx=Minus_Intersections

                ray_ind, plus_range_ind, extra_bins = return_pnt_index(Data, gate_x, gate_y, sweep, given_x_lon=Plus_Intx.x, given_y_lat=Plus_Intx.y)
                ray_ind, minus_range_ind, extra_bins = return_pnt_index(Data, gate_x, gate_y, sweep, given_x_lon=Minus_Intx.x, given_y_lat=Minus_Intx.y)
                #iterate over each range bin
                for c_rangebin in range(len(range_bins_alongray)): 
                    if c_rangebin in np.arange(minus_range_ind, plus_range_ind):
                        radarbin_entry=subset_orig_rmom[c_ray, int(c_rangebin)]
                        if radarbin_entry== -9999.0 or subset_orig_rmom._mask[c_ray, int(c_rangebin)] == True: #if the bin is masked
                            Surge_Area_Rdata[c_ray, int(c_rangebin)] = np.nan
                        else:
                            Surge_Area_Rdata[c_ray, int(c_rangebin)] = radarbin_entry 
                    else:
                        Surge_Area_Rdata[c_ray, int(c_rangebin)] = np.nan

            else:
                Surge_Area_Rdata[c_ray, :] = np.nan


            #  if c_ray in rays[:]:
            if c_ray in np.arange(np.min(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray']),
                                  np.max(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray'])):
                ray_slope, ray_yintercept, r_value, p_value, std_err = stats.linregress(gate_x[c_ray],gate_y[c_ray])
                if config.Surge_controls['Surges']['ray_selection_method']== 'average':
                    #using the average surge slope 
                    x_crsspnt, y_crsspnt = pnt_intersection(ray_yintercept, surge_yintercept, ray_slope, surge_slope)
                    pnt_along_surge= surge_linepoint
                elif config.Surge_controls['Surges']['ray_selection_method']== 'local':
                    #using the slope of the surge at the local point 
                    try: 
                        neighbor_ray_line_object=LineString(list(zip(gate_x[c_ray+1], gate_y[c_ray+1])))
                        Neighborpnt_Intersections=Surge_Holder[Surge_ID]['Center_pnts']['line_object'].intersection(neighbor_ray_line_object)
                    except:
                        neighbor_ray_line_object=LineString(list(zip(gate_x[c_ray-1], gate_y[c_ray-1])))
                        Neighborpnt_Intersections=Surge_Holder[Surge_ID]['Center_pnts']['line_object'].intersection(neighbor_ray_line_object)
                    ray_line_object=LineString(list(zip(gate_x[c_ray], gate_y[c_ray])))
                    Intersections=Surge_Holder[Surge_ID]['Center_pnts']['line_object'].intersection(ray_line_object)
                    if isinstance(Intersections, MultiPoint):
                        Intersect= Intersections[0]
                    else: 
                        Intersect= Intersections
                    x_crsspnt, y_crsspnt= Intersect.x, Intersect.y
                    pnt_along_surge= Neighborpnt_Intersections

                rayangle=angle3pt((0, 0), (x_crsspnt, y_crsspnt), pnt_along_surge)
                rayangle_90diff = rayangle-90

                if abs(rayangle_90diff) <= config.Surge_controls['Surges']['ray_selection_thresh']:
                    ray_selection[c_ray, :] = 300 
                    store_rays_wthin_range.append(c_ray)
                else:
                    ray_selection[c_ray, :] = rayangle_90diff
                rayangle_tosurge.append(rayangle_90diff)
                
            else: 
                ray_selection[c_ray, :] = np.nan 
        rayangle_tosurge= np.asarray(rayangle_tosurge)

        Rayangle_rmom= np.zeros(r_mom_ofintrest_data.shape)
        Rayangle_rmom[sweep_startidx:sweep_endidx+1] = ray_selection 
        Data['P_Radar'].rfile.add_field('Ray_angle_'+str(Surge_ID), {'data': Rayangle_rmom}, replace_existing=True) 

        Surge_subset_rmom= np.zeros(r_mom_ofintrest_data.shape)
        Surge_subset_rmom[sweep_startidx:sweep_endidx+1] = Surge_Area_Rdata
        Data['P_Radar'].rfile.add_field('Surge_subset_'+str(Surge_ID), {'data': Surge_subset_rmom}, replace_existing=True) 
        
        #* 8 * 8 *
        Conv, X, Y, Z, valid_ray_bool, DiffCon_Holder=calc_vdiff_con(config, Data, Surge_ID, Surge_Holder[Surge_ID], store_rays_wthin_range, sweep, subset_orig_rmom, np.max(offset_bins))
        Surge_Holder[Surge_ID]['Diff_Calcs']= DiffCon_Holder




        print('HHHHHHHJJJJJJJLLLLLLHHHHH')

        # * * *
        '''
        Max_diffholder=[]
        Area_of_interest=(Surge_Area_Rdata[np.min(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray']):np.max(
                          Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray']),:])
        for ray in range(np.shape(Area_of_interest)[0]):
            ray_of_interest=Area_of_interest[ray]
            ray_wo_nan= ray_of_interest[~np.isnan(ray_of_interest)]
            if len(ray_wo_nan) != 0:
                Max_diffholder.append(np.max(ray_wo_nan.ptp()))
        Surge_Holder[Surge_ID]['Offset_diff']={'Max_alongrays': Max_diffholder} 

        for count, dist in enumerate(config.Surge_controls['Surges']['offset_dist']):
            exact_diff_holder, exactplus_bin, exactminus_bin=[],[],[]
            maxseg_diff_holder, maxsegplus_bin, maxsegminus_bin=[],[],[]
            for ray in np.arange(np.min(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray']),
                                 np.max(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray'])):
                ray_line_object=LineString(list(zip(gate_x[ray], gate_y[ray])))

                Plus_Intersections=Surge_Holder[Surge_ID]['Offset_pnts'][dist]['Plus']['line_object'].intersection(ray_line_object)
                Minus_Intersections=Surge_Holder[Surge_ID]['Offset_pnts'][dist]['Minus']['line_object'].intersection(ray_line_object)
                if isinstance(Plus_Intersections, MultiPoint): Plus_Intx=Plus_Intersections[0]
                elif isinstance(Plus_Intersections, Point): Plus_Intx=Plus_Intersections
                if isinstance(Minus_Intersections, MultiPoint): Minus_Intx=Minus_Intersections[0]
                elif isinstance(Minus_Intersections, Point): Minus_Intx=Minus_Intersections

                ray_ind, plus_range_ind, _ = return_pnt_index(Data, gate_x, gate_y, sweep, given_x_lon=Plus_Intx.x, given_y_lat=Plus_Intx.y)
                ray_ind, minus_range_ind, _ = return_pnt_index(Data, gate_x, gate_y, sweep, given_x_lon=Minus_Intx.x, given_y_lat=Minus_Intx.y)
                
                ###Exact diff at each dist
                exactplus_radarbin_entry=subset_orig_rmom[ray_ind, plus_range_ind]
                exactminus_radarbin_entry=subset_orig_rmom[ray_ind, minus_range_ind]
                if exactplus_radarbin_entry== -9999.0 or subset_orig_rmom._mask[ray_ind, int(plus_range_ind)] == True: #if the bin is masked
                    pass
                elif exactminus_radarbin_entry== -9999.0 or subset_orig_rmom._mask[ray_ind, int(minus_range_ind)] == True: #if the bin is masked
                    pass
                else:
                    exactplus_bin.append(exactplus_radarbin_entry)
                    exactminus_bin.append(exactminus_radarbin_entry)
                    exact_diff_holder.append(abs(exactplus_radarbin_entry-exactminus_radarbin_entry))

                ###Max diff in each segment 
                if count == 0: 
                    #  test=subset_orig_rmom[ray_ind, :]
                    subsetray=subset_orig_rmom[ray_ind, minus_range_ind:plus_range_ind]
                    subsetray=subsetray[subsetray != -9999.0]

                    maxseg_diff_holder.append(np.max(subsetray.ptp()))
                else: 
                    currentplus_range_ind=plus_range_ind
                    currentminus_range_ind=minus_range_ind
                    pre_plus_range_ind=plus_index_storage
                    pre_minus_range_ind=minus_index_storage
                    
                    if pre_plus_range_ind > currentplus_range_ind:
                        EP,SP= pre_plus_range_ind, currentplus_range_ind
                    else: 
                        EP,SP= currentplus_range_ind,pre_plus_range_ind 
                    if pre_minus_range_ind > currentminus_range_ind:
                        EM,SM= pre_minus_range_ind, currentminus_range_ind
                    else: 
                        EM,SM= currentminus_range_ind,pre_minus_range_ind 

                    plus_subset=subset_orig_rmom[ray_ind, SP:EP]
                    minus_subset=subset_orig_rmom[ray_ind, SM:EM]
                    plus_subsetray=plus_subset[plus_subset!= -9999.0]
                    minus_subsetray=minus_subset[minus_subset != -9999.0]
                    if len(plus_subsetray) == 0 and len(minus_subsetray)==0:
                        pass
                    else: 
                        subsetray= np.ma.append(plus_subsetray, minus_subsetray)
                        maxseg_diff_holder.append(np.max(subsetray.ptp()))

                plus_index_storage=plus_range_ind
                minus_index_storage=minus_range_ind


            Surge_Holder[Surge_ID]['Offset_diff'].update({dist: exact_diff_holder}) 
            Surge_Holder[Surge_ID]['Offset_diff'].update({str(dist)+'_Max': maxseg_diff_holder}) 

        '''


        # * * *
        maxindexlist=np.where(Surge_Area_Rdata == np.nanmax(Surge_Area_Rdata))
        minindexlist=np.where(Surge_Area_Rdata == np.nanmin(Surge_Area_Rdata))
        xmax, ymax, xmin, ymin = [],[],[],[]
        for pnt in range(len(maxindexlist[0])):
            xmax.append(gate_x[maxindexlist[0][pnt], maxindexlist[1][pnt]])
            ymax.append(gate_y[maxindexlist[0][pnt], maxindexlist[1][pnt]])
        for pnt in range(len(minindexlist[0])):
            xmin.append(gate_x[minindexlist[0][pnt], minindexlist[1][pnt]])
            ymin.append(gate_y[minindexlist[0][pnt], minindexlist[1][pnt]])
        # * * *
        Surge_Holder[Surge_ID].update({'zoom_Domain': ZOOM_Dom, 
                                       'rel_ray_angle':rayangle_tosurge,
                                       'line': {'slope':surge_slope, 'intercept':surge_yintercept},
                                       'Surge_Subset': 
                                           {'R_data':Surge_Area_Rdata,
                                            'Valid_Rays': store_rays_wthin_range,
                                            'Max':{'value':np.nanmax(Surge_Area_Rdata), 'x': xmax, 'y':ymax },
                                            'Min':{'value':np.nanmin(Surge_Area_Rdata), 'x': xmin, 'y':ymin }}
                                       })

    
        ######
        Surge_stats_dict=compute_surgestats(self, day, Data, config, Surge_ID, Surge_Holder, Conv, valid_ray_bool)
        Surge_Holder[Surge_ID].update({'Stats':Surge_stats_dict})
        
        tilt_ang = Data['P_Radar'].rfile.get_elevation(sweep)
        Conv_1dplt(Z, X, Y, Surge_ID, day, Data, config, tilt_ang, valid_ray_bool, Surge_stats_dict)
        Conv_histograms(Conv, Surge_ID, day, Data, config, tilt_ang, valid_ray_bool, Surge_stats_dict)
    #  pprint.pprint(Surge_Holder)
    return Surge_Holder, gate_x, gate_y, gate_lon, gate_lat


    # find dist from point to y intercept (aka dist c)
    #  y_intercept = Point(0,Surges[Surge_ID]['line']['intercept'])
    #  triangle_side_along_surge = ray_point.distance(y_intercept)

    # a is dist from radar to y intercept
    #  triangle_side_along_NS = Surges[Surge_ID]['line']['intercept']

    # b is dist from radar to point
    #  triangle_side_along_ray = Surges[Surge_ID]['range_bins']['center'][pnt]

    # C is the angle from radar beam to true North 
    # law of cos: cos(A) = (b^2 + c^2 -a^2)/(2bc)
    #  Numerator = (triangle_side_along_surge**2) + (triangle_side_along_ray**2) - (triangle_side_along_NS**2)
    #  Denominator = 2 * triangle_side_along_ray * triangle_side_along_surge
    #  Fraction = Numerator / Denominator
    #  result_angle = math.acos(Fraction)
    #  return result_angle

#################################################
def calc_vdiff_con(config, Data, Surge_ID, Surgedata, store_rays_wthin_range, sweep, subset_orig_rmom, Max_obins):
    gate_x, gate_y, _ = Data['P_Radar'].rfile.get_gate_x_y_z(sweep)
    range_bins_alongray = Data['P_Radar'].rfile.range['data']
    
    #set up blank lists
    Rays_conv, Rays_dist, Rays_diff=[], [],[]
    Rays_mdata, Rays_pdata=[],[]
    All_conv_mean, All_diff_mean, All_conv_max, All_diff_max=[],[],[],[]
    Vray_conv_mean, Vray_diff_mean, Vray_conv_max, Vray_diff_max=[],[],[],[]
    valid_ray_bool=[]

    #  for ray in range(np.shape(Area_of_interest)[0]):
    for ray in np.arange(np.min(Surgedata['Center_pnts']['Index']['Ray']),
                         np.max(Surgedata['Center_pnts']['Index']['Ray'])):
        #Determine the radarbin of the surge for a given ray
        ray_line_object=LineString(list(zip(gate_x[ray], gate_y[ray])))
        Intersections=Surgedata['Center_pnts']['line_object'].intersection(ray_line_object)
        if isinstance(Intersections, MultiPoint): Intersect= Intersections[0]
        else:  Intersect= Intersections
        ray_ind, c_range_ind, _ = return_pnt_index(Data, gate_x, gate_y, sweep, given_x_lon=Intersect.x, given_y_lat=Intersect.y)
        exact_c_radarbin_entry=subset_orig_rmom[ray_ind, c_range_ind]

        m_rangebound=c_range_ind-Max_obins
        p_rangebound=c_range_ind+Max_obins
        if m_rangebound <= 0:
            Mdata=subset_orig_rmom[ray_ind, 0:c_range_ind]
            while len(Mdata) < Max_obins:
                Mdata=np.append(Mdata,np.nan)
            #  Mdata=ma.full(Max_obins, np.nan)
            #  Mdata[0:Max_obins+1]=subset_orig_rmom[ray_ind, 0:c_range_ind]
        else:
            Mdata=subset_orig_rmom[ray_ind, m_rangebound:c_range_ind]

        Mdata[Mdata==-9999.0]=np.nan
        Mdata[Mdata._mask==True]=np.nan
        Pdata=subset_orig_rmom[ray_ind, c_range_ind:p_rangebound]
        Pdata[Pdata==-9999.0]=np.nan
        Pdata[Pdata._mask==True]=np.nan

        Convergence, Difference, Distance=[], [],[]
        for obin in range(1, Max_obins+1):
            dist= (obin*Data['P_Radar'].rfile.range['meters_between_gates'])*2
            Distance.append(dist)
            plus_range_ind, minus_range_ind = (c_range_ind + obin), (c_range_ind - obin)

            if plus_range_ind <= len(range_bins_alongray):
                p_rb_entry, p_rb_entry_m = subset_orig_rmom[ray_ind, plus_range_ind], subset_orig_rmom._mask[ray_ind, plus_range_ind]
                m_rb_entry, m_rb_entry_m = subset_orig_rmom[ray_ind, minus_range_ind], subset_orig_rmom._mask[ray_ind, minus_range_ind]

                if p_rb_entry== -9999.0 or p_rb_entry_m == True or m_rb_entry== -9999.0 or m_rb_entry_m == True: #if either bin is masked
                    Convergence.append(np.nan)
                    Difference.append(np.nan)
                else:#actually calc convergence 
                    diff = m_rb_entry - p_rb_entry
                    Difference.append(diff)
                    con = diff/dist
                    Convergence.append(con)

        All_conv_max.append(np.nanmax(Convergence))
        All_diff_max.append(np.nanmax(Difference))

        Count=np.count_nonzero(~np.isnan(Difference))
        Ave_Conv=np.nanmean(Convergence)
        Ave_Diff=np.nanmean(Difference)
        All_conv_mean.append(Ave_Conv)
        All_diff_mean.append(Ave_Diff)
        #  print(Ave_Conv)
        #  print('AveCon: {}\nCount: {}\n\n'.format(Ave_Conv, Count))
        #  print('Con: {}\nDiff: {}\nDist: {}\nCount: {}\n\n'.format(Convergence, Difference, Distance, Count))
        Rays_pdata.append(Pdata)
        Rays_mdata.append(Mdata)
        Rays_conv.append(Convergence)
        Rays_dist.append(Distance)
        Rays_diff.append(Difference)

        if ray in store_rays_wthin_range:
            Vray_conv_max.append(np.nanmax(Convergence))
            Vray_diff_max.append(np.nanmax(Difference))
            Vray_conv_mean.append(Ave_Conv)
            Vray_diff_mean.append(Ave_Diff)
            valid_ray_bool.append(True)
        else:
            valid_ray_bool.append(False)
        
    # * * * *
    Conv, Dist, Diff, Md, Pd = np.array(Rays_conv), np.array(Rays_dist), np.array(Rays_diff), np.array(Rays_mdata), np.array(Rays_pdata)
    Z=[Conv, Diff, Md, Pd, Dist ]
    X,Y=np.meshgrid(range(Conv.shape[1]+1), range(Conv.shape[0]+1))
    # * * * *

    DiffCon_Holder={'Conv':{'Max':{'Vrays':Vray_conv_max, 'Allrays':All_conv_max},
                            'Mean':{'Vrays':Vray_conv_mean, 'Allrays':All_conv_mean}, 
                            'PTP':{'Vrays':None, 'Allrays':None}},
                    'Vdiff_convcalc':{'Max':{'Vrays':Vray_diff_max, 'Allrays':All_diff_max},
                                      'Mean':{'Vrays':Vray_diff_mean, 'Allrays':All_diff_mean}, 
                                      'PTP':{'Vrays':None, 'Allrays':None}},
                    'Vdiff_odist':{'Max':{'Vrays':None, 'Allrays':None},
                                   'Mean':{'Vrays':None, 'Allrays':None}, 
                                   'PTP':{'Vrays':None, 'Allrays':None}}
                    }

    return Conv, X,Y,Z, valid_ray_bool, DiffCon_Holder




