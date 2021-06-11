#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import needed modules
######################
import matplotlib
#  matplotlib.use('agg')
#  matplotlib.use('GTK3Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv 
import math 
import pprint
from scipy import stats
from scipy import ndimage, interpolate
from shapely.geometry import LineString, Point, Polygon, MultiPoint

#this is the file with the plotting controls to access any of the vars in that file use config.var
import config as plot_config
from read_pforms import Platform

def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

def clicker_defn(PLT):
    pts = [np.nan]
    #  plt.setp(plt.gca(), autoscale_on=False)
    #  plt.setp(plt.gcf())
    #  plt.ion()
    #  plt.draw()
    # Now do a zoom
    tellme('Now do a nested zoom, click to begin')
    plt.waitforbuttonpress()

    while True:
        tellme('Select two corners of zoom, middle mouse button to finish')
        pts = plt.ginput(2, timeout=-1)
        if len(pts) < 2:
            break
        (x0, y0), (x1, y1) = pts
        xmin, xmax = sorted([x0, x1])
        ymin, ymax = sorted([y0, y1])
        AXES=PLT.fig.get_axes()
        for a in AXES[::2]: 
                plt.sca(a)
                plt.xlim(xmin, xmax)
                plt.ylim(ymin, ymax)
    tellme('I SWEAR YOU BETTER WORK THIS TIME....  OR THIS RELATIONSHIP IS OVER!!!')
    #  plt.waitforbuttonpress()
    pts = plt.ginput(n=40,timeout=0) # look up ginput docs for some good guidance
    #  plt.waitforbuttonpress()
    return pts

def make_csv(Data, pts, surge_name, filename, Does_csv_exist):
    if Does_csv_exist == False: 
        # making csv file if it does not exist 
        with open(filename, 'w') as csvfile: 
            # field names 
            fields = ['RName','Sweep', 'Scan_time','Surge_ID', 'point_x', 'point_y'] 
            
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

            rows= [[Data['P_Radar'].name, Data['P_Radar'].swp, Data['P_Radar'].Scan_time, surge_name, point[0],point[1]]]
            csvwriter.writerows(rows)

def return_pnt_index(self, Data, surge_pd, gate_x, gate_y, sweep, Buffer_dist=None, point=None, given_x=None, given_y=None):
    def find_index_of_nearest_xy(y_array, x_array, y_point, x_point):
        distance = (y_array-y_point)**2 + (x_array-x_point)**2
        idy,idx = np.where(distance==distance.min())
        return idy[0],idx[0]
    # * * * *
    if point != None:
        pnt_X, pnt_Y= surge_pd.point_x[point], surge_pd.point_y[point]
    else: 
        pnt_X, pnt_Y = given_x, given_y 
    ray_ind, range_ind= find_index_of_nearest_xy(gate_y, gate_x, pnt_Y, pnt_X)

    if Buffer_dist != None:
        bin_dist= Data['P_Radar'].rfile.range['meters_between_gates']
        extra_bins = round(Buffer_dist/bin_dist)
    else: extra_bins = None
    #  self.display.plot_line_xy(surge_pd.point_x, surge_pd.point_y)
    #  pnt=points.split("), (")
    return ray_ind, range_ind, extra_bins



def surge_csv_reader(config, Data, day):        
    def subset_surge_df(surge_df, Data, ID):
        surge_ID_sub=surge_df[surge_df['Surge_ID'] == ID]
        time_sub=surge_ID_sub[surge_ID_sub['Scan_time'] == str(Data['P_Radar'].Scan_time)]
        radar_sub=time_sub[time_sub['RName']==Data['P_Radar'].name]
        final_sub=radar_sub[radar_sub['Sweep']==Data['P_Radar'].swp]
        final_sub.reset_index()
        return final_sub
    # ****

    surge_df=pd.read_csv(config.g_TORUS_directory+day+'/data/'+day+'_surge_pnts.csv')
    surges=surge_df.Surge_ID.unique()

    num_of_surges, surge_ids, surges_attime_df = 0, [], pd.DataFrame() 
    for j in surges: 
        single_surge_df=subset_surge_df(surge_df, Data, j)
        if len(single_surge_df.point_x) == 0:
            pass
        else:
            surges_attime_df=surges_attime_df.append(single_surge_df)
            num_of_surges=num_of_surges+1
            surge_ids.append(j)
    return num_of_surges, surge_ids, surges_attime_df  


def surge_radarbins(self, config, Data, sweep):
    def angle3pt(a, b, c):
        """Counterclockwise angle in degrees by turning from a to c around b
            Returns a float between 0.0 and 360.0"""
        ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
        if b in (a, c):
            raise ValueError("Undefined angle, two identical points", (a, b, c))
        return ang + 360 if ang < 0 else ang
    def pnt_intersection(b1, b2, m1, m2):
        xi = (b1-b2) / (m2-m1)
        yi = m1*xi + b1
        return xi, yi
    
    gate_x, gate_y, _ = Data['P_Radar'].rfile.get_gate_x_y_z(sweep)
    gate_lat, gate_lon, _ = Data['P_Radar'].rfile.get_gate_lat_lon_alt(sweep)
    range_bins_alongray = Data['P_Radar'].rfile.range['data']
    r_mom_ofintrest_data=Data['P_Radar'].rfile.fields[config.lineplt_control['H']['var']]['data']
    sweep_startidx = np.int64(Data['P_Radar'].rfile.sweep_start_ray_index['data'][sweep])
    sweep_endidx = np.int64(Data['P_Radar'].rfile.sweep_end_ray_index['data'][sweep])

    Surge_Holder={}
    for Surge_ID in self.surge_ids[:]:
        pnt_name = []
        rays_ind, ranges_ind = [], []
        r_bin_sub, r_data_sub = [], []

        zoom_surge_df= self.valid_surges[self.valid_surges['Surge_ID'] == Surge_ID]
        for i in zoom_surge_df.index.tolist():
            pnt_name.append(i)
            #determine the index positions for each point along the surge (and for each offset point)
            offset_bins=[]
            for offset_dist in config.Surge_controls['offset_dist']:
                ray_ind, range_ind, extra_bins = return_pnt_index(self, Data, self.valid_surges, gate_x, gate_y, sweep, Buffer_dist=offset_dist, point=i)
                offset_bins.append(extra_bins)
            rays_ind.append(ray_ind)
            ranges_ind.append(range_ind)
        Surge_Holder.update({Surge_ID:{'point_labels': pnt_name}})  
        Surge_Holder[Surge_ID].update({'Offset_pnts': 
                                            {'Max': {'offset': np.max(offset_bins),'dist': np.max(config.Surge_controls['offset_dist'])},
                                             'Min': {'offset': np.min(offset_bins),'dist': np.min(config.Surge_controls['offset_dist'])}
                                             }
                                        })

           
        #store the relevant data for each surge (and surge offset) point in a dict 
        for o in range(len(config.Surge_controls['offset_dist'])):
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
                    Surge_Holder[Surge_ID]['Offset_pnts'].update({config.Surge_controls['offset_dist'][o]: Both_Dir}) 
                
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

        #add Polygon for the full (maxoffset) area 
        Max_odist=np.max(config.Surge_controls['offset_dist'])
        Xplus=Surge_Holder[Surge_ID]['Offset_pnts'][Max_odist]['Plus']['x']
        Yplus=Surge_Holder[Surge_ID]['Offset_pnts'][Max_odist]['Plus']['y']
        Xminus=Surge_Holder[Surge_ID]['Offset_pnts'][Max_odist]['Minus']['x']
        Yminus=Surge_Holder[Surge_ID]['Offset_pnts'][Max_odist]['Minus']['y']
        
        comb_xy_minus, comb_xy_plus = list(zip(Xminus, Yminus)), list(zip(Xplus, Yplus))
        poly=Polygon(comb_xy_plus+comb_xy_minus)
        Surge_Holder[Surge_ID].update({'Polygon': poly}) 

        # * * *
        if 'zoom' in config.r_mom: 
            c_x, c_y = zoom_surge_df.point_x.mean(), zoom_surge_df.point_y.mean()
            cray_ind, crange_ind, _ = return_pnt_index(self, Data, zoom_surge_df, gate_x, gate_y, sweep, given_x=c_x, given_y=c_y)
            p_lat, p_lon = gate_lat[cray_ind, crange_ind], gate_lon[cray_ind, crange_ind]

            if self.config.radar_controls['offsetkm']['Zoom'] == 'Fitted':
                surge_extents=Surge_Holder[Surge_ID]['Polygon'].envelope.bounds
                x_dist=(surge_extents[2]-surge_extents[0])/2000
                y_dist=(surge_extents[3]-surge_extents[1])/2000
                if x_dist > y_dist: o_km=x_dist+.4 #.1 is to add an additional buffer
                elif y_dist >= x_dist: o_km=y_dist+.4
            else:
                o_km=self.config.radar_controls['offsetkm']['Zoom']

            ZOOM_Dom = Platform.getLocation('ClickPoint', scan_time=Data['P_Radar'].Scan_time, 
                                              offsetkm=o_km , pnt_x=p_lon, pnt_y=p_lat)
        else: 
            ZOOM_Dom = None

        # * * *
        points = [Surge_Holder[Surge_ID]['Center_pnts']['line_object'].interpolate(dist) for dist in 
                  np.arange(0, Surge_Holder[Surge_ID]['Center_pnts']['line_object'].length, 0.5)]
        x, y = zip(*[(pt.x, pt.y) for pt in points])
        surge_slope, surge_yintercept, r_value, p_value, std_err = stats.linregress(x,y)

        x_vals = (np.min(gate_x), np.max(gate_x))
        x_vals= np.asarray(x_vals)
        y_vals = surge_yintercept + surge_slope * x_vals
        surge_linepoint= (x_vals[-1], y_vals[-1])

        # * * *
        #define the indices for the required sweep
        ray_selection= np.zeros(gate_x.shape)
        rayangle_tosurge=[]
        subset_orig_rmom=r_mom_ofintrest_data[sweep_startidx:sweep_endidx+1]
        #  Surge_Area_Rdata=r_mom_ofintrest_data[sweep_startidx:sweep_endidx+1]
        Surge_Area_Rdata=np.zeros(gate_x.shape)
        for c_ray in range(ray_selection.shape[0]): #iterate over each ray
            if c_ray in np.arange(np.min(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray']),
                                  np.max(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray'])):
                ray_line_object=LineString(list(zip(gate_x[c_ray], gate_y[c_ray])))
                Max_dist = Surge_Holder[Surge_ID]['Offset_pnts']['Max']['dist']
                Plus_Intersections=Surge_Holder[Surge_ID]['Offset_pnts'][Max_dist]['Plus']['line_object'].intersection(ray_line_object)
                Minus_Intersections=Surge_Holder[Surge_ID]['Offset_pnts'][Max_dist]['Minus']['line_object'].intersection(ray_line_object)
                Plus_ind, Minus_ind =[],[]

                if isinstance(Plus_Intersections, Point):
                    #  print(type(Plus_Intersections))
                    ray_ind, plus_range_ind, extra_bins = return_pnt_index(self, Data, self.valid_surges, gate_x, gate_y, sweep, given_x=Plus_Intersections.x, given_y=Plus_Intersections.y)
                    ray_ind, minus_range_ind, extra_bins = return_pnt_index(self, Data, self.valid_surges, gate_x, gate_y, sweep, given_x=Minus_Intersections.x, given_y=Minus_Intersections.y)
                    #iterate over each range bin
                    for c_rangebin in range(len(range_bins_alongray)): 
                        if Surge_Area_Rdata[c_ray, int(c_rangebin)] == -9999.0:
                            Surge_Area_Rdata[c_ray, int(c_rangebin)] = np.nan
                        
                        if c_rangebin in np.arange(minus_range_ind, plus_range_ind):
                            Surge_Area_Rdata[c_ray, int(c_rangebin)] = subset_orig_rmom[c_ray, int(c_rangebin)] 
                        else:
                            #  Surge_Area_Rdata[c_ray, int(c_rangebin)] = 0
                            Surge_Area_Rdata[c_ray, int(c_rangebin)] = np.nan
            else:
                Surge_Area_Rdata[c_ray, :] = np.nan


            #  if c_ray in rays[:]:
            if c_ray in np.arange(np.min(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray']),
                                  np.max(Surge_Holder[Surge_ID]['Center_pnts']['Index']['Ray'])):
                ray_slope, ray_yintercept, r_value, p_value, std_err = stats.linregress(gate_x[c_ray],gate_y[c_ray])
                xi, yi = pnt_intersection(ray_yintercept, surge_yintercept, ray_slope, surge_slope)
                rayangle=angle3pt((0, 0), (xi, yi), surge_linepoint)
                rayangle_90diff = rayangle-90

                if abs(rayangle_90diff) <= config.Surge_controls['ray_selection']:
                    ray_selection[c_ray, :] = 300 
                else:
                    ray_selection[c_ray, :] = rayangle_90diff
                rayangle_tosurge.append(rayangle_90diff)
            else: 
                ray_selection[c_ray, :] = np.nan 
        #  Surge_Area_Rdata = np.asarray(Surge_Area_Rdata)
        rayangle_tosurge= np.asarray(rayangle_tosurge)

        Rayangle_rmom= np.zeros(r_mom_ofintrest_data.shape)
        Rayangle_rmom[sweep_startidx:sweep_endidx+1] = ray_selection 
        Data['P_Radar'].rfile.add_field('Ray_angle', {'data': Rayangle_rmom}, replace_existing=True) 

        #  print(np.shape(Surge_Area_Rdata))
        Surge_subset_rmom= np.zeros(r_mom_ofintrest_data.shape)
        #  Surge_subset_rmom= Surge_Area_Rdata
        Surge_subset_rmom[sweep_startidx:sweep_endidx+1] = Surge_Area_Rdata
        Data['P_Radar'].rfile.add_field('Surge_subset_'+str(Surge_ID), {'data': Surge_subset_rmom}, replace_existing=True) 

        #  print(np.where(x==np.min(Surge_Area_Rdata)))
        #  print(np.amin(Surge_Area_Rdata))
        #  result = np.where(Surge_subset_rmom== np.nanmin(Surge_subset_rmom))
        result= np.nanmin(Surge_subset_rmom)
        #  print(result[0])
        #  print(result[1])
        #  test=gate_x(result[0], result[1])
        #  print(test)
        # * * *
        Surge_Holder[Surge_ID].update({'zoom_Domain': ZOOM_Dom, 
                                       'rel_ray_angle':rayangle_tosurge,
                                       'line': {'slope':surge_slope, 'intercept':surge_yintercept},
                                       'Surge_Subset': 
                                           {'R_data:':Surge_Area_Rdata, 
                                            'Max':{'value':np.nanmax(Surge_Area_Rdata), 'x': None, 'y':None },
                                            'Min':{'value':np.nanmin(Surge_Area_Rdata), 'x': None , 'y':None }}
                                       })
    #  pprint.pprint(Surge_Holder)
    print('made it into here')
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

##################################################
# Define a triangle by clicking three points
'''
plt.clf()
plt.setp(plt.gca(), autoscale_on=False)

tellme('You will define a triangle, click to begin')
plt.waitforbuttonpress()

while True:
    pts = []
    while len(pts) < 3:
        tellme('Select 3 corners with mouse')
        pts = np.asarray(plt.ginput(3, timeout=-1))
        if len(pts) < 3:
            tellme('Too few points, starting over')
            time.sleep(1)  # Wait a second

    ph = plt.fill(pts[:, 0], pts[:, 1], 'r', lw=2)

    tellme('Happy? Key click for yes, mouse click for no')

    if plt.waitforbuttonpress():
        break

    # Get rid of fill
    for p in ph:
        p.remove()


##################################################
# Now contour according to distance from triangle
# corners - just an example

# Define a nice function of distance from individual pts
def f(x, y, pts):
    z = np.zeros_like(x)
    for p in pts:
        z = z + 1/(np.sqrt((x - p[0])**2 + (y - p[1])**2))
    return 1/z


X, Y = np.meshgrid(np.linspace(-1, 1, 51), np.linspace(-1, 1, 51))
Z = f(X, Y, pts)

CS = plt.contour(X, Y, Z, 20)

tellme('Use mouse to select contour label locations, middle button to finish')
CL = plt.clabel(CS, manual=True)

##################################################
# Now do a zoom
tellme('Now do a nested zoom, click to begin')
plt.waitforbuttonpress()

while True:
    tellme('Select two corners of zoom, middle mouse button to finish')
    pts = plt.ginput(2, timeout=-1)
    if len(pts) < 2:
        break
    (x0, y0), (x1, y1) = pts
    xmin, xmax = sorted([x0, x1])
    ymin, ymax = sorted([y0, y1])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

tellme('All Done!')
plt.show()

'''

