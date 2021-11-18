#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse, cProfile, logging, time, os, os.path
import pprint
import functools
import hdbscan
import mpu
import time
import json
import pickle
import pprint 
import jsonpickle
import jsonpickle.ext.numpy as jsonpickle_numpy
from functools import wraps
from scipy import spatial
from pathlib import Path
from tabulate import tabulate
from sklearn.cluster import DBSCAN
from itertools import combinations, permutations

from read_pforms import read_TInsitu, pform_names
import config as plot_config

jsonpickle_numpy.register_handlers()

#pdt=lambda df:tabulate(df,headers='keys',tablefmt='psql')
pdt=lambda df:tabulate(df, headers='keys',tablefmt='pretty', floatfmt=".4f")

#width = pd.util.terminal.get_terminal_size() # find the width of the user's terminal window
#pd.set_option('display.width', width[0]) # set that as the max width in Pandas
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.options.display.width=None
np.set_printoptions(precision=6)


## Thresholds
# # # # # # #
#frequency at which to bin (average) the obs
#  bined_data_freq = '10S'
#threshold between two obs for valid comparisions
#  valid_time_diff= 5 #Min
valid_dist_diff= 10 #'10km'
# parameters to cross compare
parameters_to_difference = ['Thetav', 'Thetae', 'tfast']
stats_to_compute=['mean','std','count']

day_list= ['20190517','20190518','20190520', '20190523', '20190524','20190525','20190526',
             '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
#  day_list= ['20190517','20190518','20190520']
#  day_list= ['20190517','20190518','20190520', '20190524','20190525','20190526',
             #  '20190527', '20190528', '20190608', '20190611', '20190613', '20190615']
day_cuttoff= {'20190517':'23:00', '20190518':'22:30','20190520':'17:00', '20190523':'22:00', '20190524':'19:00','20190525':'19:00','20190526':'18:00',
              '20190527':'20:00', '20190528':'22:15', '20190608':'20:30', '20190611':'22:50', '20190613':'23:59', '20190615':'21:30'}#'22:00'}
day_startoff = {'20190615': '20:30'}
colors ={'Prb1': 'red', 'Prb2':'green', 'WinS':'blue', 'LIDR': 'black', 'FFld':'gray', 'CoMeT1':'purple', 'CoMeT2':'pink', 'CoMeT3':'orange'}

perms=list(permutations(pform_names('TInsitu'), 2))
col_names=[]
for c in range(len(perms)):
    #  print('Perms: {}, {} NCS:{}'.format(perms[c][0], perms[c][1], NCS))
    col_names.append(perms[c][0]+'-'+perms[c][1])
iterables=[day_list, parameters_to_difference, stats_to_compute]
Multiday_df_index= pd.MultiIndex.from_product(iterables, names=["day", "parameter", "stat"])
Multiday_df=pd.DataFrame(data=np.nan, index=Multiday_df_index, columns=col_names)

#######################################################

def timeit(my_func):
    @wraps(my_func)
    def timed(*args, **kw):

        tstart = time.time()
        output = my_func(*args, **kw)
        tend = time.time()

        print('"{}" took {:.3f} ms to execute\n'.format(my_func.__name__, (tend - tstart) * 1000))
        return output
    return timed

@timeit
def gather_measurements(day, plot_config):
    measurement_sets = []

    platform_names = pform_names('TInsitu')


    for platform_name in platform_names:

        d_avail=read_TInsitu(plot_config, day, platform_name, True)

        if d_avail == False:
            continue

        print("processing day:", day, " platform:", platform_name)

        df, platform_type = read_TInsitu(plot_config, day, platform_name)

        #remove unneeded data columns
        df = df.drop(columns=['U','V', 'dir', 'spd' ])

        if platform_type =='NSSL':
            df = df[df['qc1']==0]
            df = df[df['qc4']==0]
            
            #convert to K
            df.tfast=df.tfast+273.15

            # More agressive quality controll
            # df = df[df['all_qc_flags']==0]

            df = df.drop(columns=['id', 'qc1', 'qc2', 'qc3', 'qc4', 'all_qc_flags'])

        #drop rows with no data
        df.dropna(subset=['lat', 'lon'], inplace=True)

        print(df)
        #add column with the instrument name
        df.loc[:,'pname'] = platform_name
        df.loc[:,'ptype'] = platform_type

        measurement_sets.append(df)

    # Put all measuremens from all devices on all days in one big measurement data frame
    measurements = pd.concat(measurement_sets)

    measurements = measurements.set_index('datetime')
    measurements = measurements.sort_index()

    print("len measurements", len(measurements))

    return measurements, platform_names

@timeit
def file_to_measurements(filename = "time_sorted_measurements.json"):
    print("reading json", filename)
    with open(filename, 'r') as f:
        data = f.read()
        measurements = jsonpickle.decode(data)
    return measurements


# Find lat/lon central to all instruments in dataframe
# Add distance from instrument to center point for each instrument
def add_distance_to_center_point(df):

    lat_lon = df[['lat','lon']].to_numpy()
    center_lat, center_lon = df['lat'].mean(), df['lon'].mean()

    dist = []
    for lat,lon in lat_lon:
        dist.append( mpu.haversine_distance((center_lat, center_lon), (lat, lon)))

    df['km_to_center'] = dist

    return df


def cluster_locations_in_sample(sample, radius_km):

    #clusterer = hdbscan.HDBSCAN( cluster_selection_epsilon=radius_km/6371.0, min_cluster_size=2, metric='haversine')
    clusterer = DBSCAN(eps=radius_km/6371.0, min_samples=2)

    lat_lon = sample[['lat','lon']].to_numpy()
    rlat_lon = np.radians(lat_lon)
    try:
        clusterer.fit(rlat_lon)
        sample_set = sample.assign(cluster= clusterer.labels_)
    except:
        sample_set = sample.assign(cluster= -1)

    return sample_set


def split_clusters(sample):
    clustered_samples = []

    grouped = sample.groupby('cluster')
    for label, group in grouped:

        if label < 0:
            continue    # Skip Unclustered Instruments

        clustered_samples.append(group)

    return clustered_samples


def assign_unique_location_id(time_location_groups):
    location_id = 0

    for sample in time_location_groups:
        sample['location_id'] = location_id
        location_id += 1

    return time_location_groups


@timeit
def approach_single_data_frame (day, plot_config, outdir, parameters_to_difference, Multiday_df):

    # Gather all instruments into one big dataframe
    measurements, platform_names = gather_measurements(day, plot_config)

    # Option to dump to file / read to file if above takes a long time
    #measurements_to_file(measurements, filename = "time_sorted_measurements.json")
    #measurements = file_to_measurements(filename = "time_sorted_measurements.json")

    # Subset time during testing
    #measurements = measurements.between_time('23:38', '23:42')
    #only will consider this time range #*************
    # ******************!!!!!!!!!!!!!!!!
    if day =='20190615':
        start=day_startoff[day]
        print('day_startoff')
    else:
        start='00:01'
    startsplit=start.split(':')
    endsplit=day_cuttoff[day].split(':')
    timerange_str=startsplit[0]+startsplit[1]+'_'+endsplit[0]+endsplit[1]

    measurements = measurements.between_time(start, day_cuttoff[day])

    #
    # Resample all samples to one minute intervals.
    # This greatly reduces the data size
    # Instruments are only included if they sampled during the interval
    # Each instrument is provides 0 or 1 sample every period
    #
    instruments =  measurements.groupby('pname').resample('1min').median()
    instruments = instruments.sort_index(level='datetime')

    with open(outdir+'instr_'+day+'_'+timerange_str+'.txt', 'w') as f:
        print("instruments\n################################", file=f)
        print(pdt(instruments), file=f)
        print("\n\ninstruments.index\n################################",file=f)
        print(instruments.index, file=f)
    f.close()

    # instruments now contains a much smaller set of samples.
    # each sample is for a single instrument
    # key here...  we resampled to minute granularity
    # so if you groupby (below) datetime you get a group of unique instruments at a one minue time stamp.
    # Dump the data structure to be clear.
        
    time_location_groups = []

    with open(outdir+'instr_timeslice_'+day+'_'+timerange_str+'.txt', 'w') as f:
        for time, instruments_at_time in instruments.groupby(level='datetime'):

            # assign cluster id to each instrument measurement
            # ******************!!!!!!!!!!!!!!!!
            radius_km = valid_dist_diff
            instruments_at_time = cluster_locations_in_sample(instruments_at_time, radius_km)
            print("\n\ninstruments at time slice", file=f)
            print(instruments_at_time, file=f)

            # split instruments into groups based on cluster id
            # location_groups:
            # [
            #   [instruments clustered at location 1],
            #   [instruments clustered at location 2],
            #   [instruments clustered at location n]
            # ]
            location_groups = split_clusters(instruments_at_time)

            for sample in location_groups:
                time_location_groups.append(sample)
    f.close()


    time_location_groups = assign_unique_location_id(time_location_groups)
    # time_location_groups = [
    #   [instruments clustered at location 1, time 1],
    #   [instruments clustered at location 2, time 1],
    #   [instruments clustered at location 3, time 1],
    #   [instruments clustered at location 4, time 2],
    #   [instruments clustered at location 5, time 2],
    #   [instruments clustered at location 6, time 3],
    # ]


    # Add column for distance to location id center point for sanity checking
    for sample in time_location_groups:
        sample = add_distance_to_center_point(sample)

    # Measurement Data Frame
    mdf = pd.concat(time_location_groups)
    #
    # Each data frame row is an instrument measurement
    #   time is reduced to samples every x seconds
    #   location is reduced to location_id's where all measurements at a location_id are within x km appart
    #
    # mdf = [
    #     [time = 0, location_id = 0, km_to_location_center, instrument = i1, lat, lon, temp],
    #     [time = 0, location_id = 0, km_to_location_center, instrument = i2, lat, lon, temp],
    #
    #     [time = 0, location_id = 1, km_to_location_center, instrument = i3, lat, lon, temp],
    #     [time = 0, location_id = 1, km_to_location_center, instrument = i4, lat, lon, temp],
    #
    #     [time = 1, location_id = 2, km_to_location_center, instrument = i1, lat, lon, temp],
    #     [time = 1, location_id = 2, km_to_location_center, instrument = i4, lat, lon, temp],
    #    ...
    # ]

    #  print("\nMeasurement Data Frame", file=f)
    #  print(mdf, file=f)

    # At this point we have a nice reduced set of times and locations
    # Do the bias calcs
    #

    new_groups = []
    for location_id, measurements_at_location in mdf.groupby('location_id'):

        instruments_at_location = measurements_at_location.unstack().index

        #  print(location_id)
        #  print(instruments_at_location)
        #  print(measurements_at_location)
        for parameter in parameters_to_difference:

            # Matrix of differences between parameters:
            # [
            #             p0:       p1:     p2:
            #     p0: [ p0 - p0, p1 - p0, p2 - p0],
            #     p1: [ p0 - p1, p1 - p1, p2 - p1],
            #     p2: [ p0 - p2, p1 - p2, p2 - p2]
            # ]
            diff_matrix = measurements_at_location[parameter].values - measurements_at_location[parameter].values[:, None]


            iterables=[[parameter], instruments_at_location]
            new_columns = pd.MultiIndex.from_product(iterables, names=["var", "delta"])
            ddf = pd.DataFrame(index=measurements_at_location.index, columns=new_columns, data=diff_matrix)

            # ddf(delta dataframe) looks like this:
            #
            # var                           Thetav
            # delta                         CoMeT1    CoMeT2    CoMeT3
            # pname  datetime
            # CoMeT1 2019-06-15 20:15:00  0.000000 -0.320648  0.071579
            # CoMeT2 2019-06-15 20:15:00  0.320648  0.000000  0.392227
            # CoMeT3 2019-06-15 20:15:00 -0.071579 -0.392227  0.000000


            measurements_at_location = measurements_at_location.join(ddf)
        new_groups.append(measurements_at_location)

    mdf = pd.concat(new_groups)

    # Master data frame rows now report the difference of "params" between other instruments at the same time and location
    #
    #                                   lat         lon       tfast      Thetav      Thetae  cluster  location_id  km_to_center  (Thetav, CoMeT1)  (Thetav, CoMeT2)  (Thetav, CoMeT3)  (Thetae, CoMeT1)  (Thetae, CoMeT2)  (Thetae, CoMeT3)  (tfast, CoMeT1)  (tfast, CoMeT2)  (tfast, CoMeT3)  (Thetav, FFld)  (Thetav, LIDR)  (Thetav, Prb1)  (Thetae, FFld)  (Thetae, LIDR)  (Thetae, Prb1)  (tfast, FFld)  (tfast, LIDR)  (tfast, Prb1)  (Thetav, Prb2)  (Thetae, Prb2)  (tfast, Prb2)
    # pname  datetime
    # CoMeT1 2019-06-15 20:15:00  35.184559 -101.934631  302.315002  315.717010  355.531998        0            0      0.000879          0.000000         -0.320648          0.071579          0.000000          0.144440          0.566208         0.000000        -0.384995         0.000000             NaN             NaN             NaN             NaN             NaN             NaN            NaN            NaN            NaN             NaN             NaN            NaN
    # CoMeT2 2019-06-15 20:15:00  35.184540 -101.934593  301.930008  315.396362  355.676437        0            0      0.003477          0.320648          0.000000          0.392227         -0.144440          0.000000          0.421768         0.384995         0.000000         0.384995             NaN             NaN             NaN             NaN             NaN             NaN            NaN            NaN            NaN             NaN             NaN            NaN
    # CoMeT3 2019-06-15 20:15:00  35.184555 -101.934662  302.315002  315.788589  356.098206        0            0      0.003034         -0.071579         -0.392227          0.000000         -0.566208         -0.421768          0.000000         0.000000        -0.384995         0.000000             NaN             NaN             NaN             NaN             NaN             NaN            NaN            NaN            NaN             NaN             NaN            NaN
    # CoMeT1 2019-06-15 20:16:00  35.184566 -101.934639  301.705002  315.061493  353.876495        0            1      0.001623          0.000000         -0.194290         -0.142975          0.000000          0.441589          0.541550         0.000000        -0.225006        -0.179993             NaN             NaN             NaN             NaN             NaN             NaN            NaN            NaN            NaN             NaN             NaN            NaN
    # CoMeT2 2019-06-15 20:16:00  35.184540 -101.934593  301.479996  314.867203  354.318085        0            1      0.003959          0.194290          0.000000          0.051315         -0.441589          0.000000          0.099960         0.225006         0.000000         0.045013             NaN             NaN             NaN             NaN             NaN             NaN            NaN            NaN            NaN             NaN             NaN            NaN
    # ...                               ...         ...         ...         ...         ...      ...          ...           ...               ...               ...               ...               ...               ...               ...              ...              ...              ...             ...             ...             ...             ...             ...             ...            ...            ...            ...             ...             ...            ...
    # CoMeT1 2019-06-15 21:15:00  35.170799 -101.938568  301.929993  315.605988  357.519012        0           60      5.358886          0.000000          0.952362               NaN          0.000000          1.998993               NaN         0.000000         0.850006              NaN        0.142574       -0.486051       -0.772723        0.157493        1.575605       -1.538068    -273.339993    -274.019993    -274.089993             NaN             NaN            NaN
    # CoMeT2 2019-06-15 21:15:00  35.170761 -101.938530  302.779999  316.558350  359.518005        0           60      5.363296         -0.952362          0.000000               NaN         -1.998993          0.000000               NaN        -0.850006         0.000000              NaN       -0.809788       -1.438413       -1.725085       -1.841500       -0.423388       -3.537061    -274.189999    -274.869999    -274.939999             NaN             NaN            NaN
    # FFld   2019-06-15 21:15:00  35.191100 -102.053900   28.590000  315.748561  357.676506        0           60      5.370685         -0.142574          0.809788               NaN         -0.157493          1.841500               NaN       273.339993       274.189999              NaN        0.000000       -0.628624       -0.915297        0.000000        1.418111       -1.695561       0.000000      -0.680000      -0.750000             NaN             NaN            NaN
    # LIDR   2019-06-15 21:15:00  35.189700 -101.992100   27.910000  315.119937  359.094617        0           60      0.843074          0.486051          1.438413               NaN         -1.575605          0.423388               NaN       274.019993       274.869999              NaN        0.628624        0.000000       -0.286673       -1.418111        0.000000       -3.113673       0.680000       0.000000      -0.070000             NaN             NaN            NaN
    # Prb1   2019-06-15 21:15:00  35.191200 -102.055400   27.840000  314.833264  355.980944        0           60      5.506871          0.772723          1.725085               NaN          1.538068          3.537061               NaN       274.089993       274.939999              NaN        0.915297        0.286673        0.000000        1.695561        3.113673        0.000000       0.750000       0.070000       0.000000             NaN             NaN            NaN

     


    with open(outdir+'mdf_'+day+'_'+timerange_str+'.txt', 'w') as f:
        print("\nMaster Data Frame", file=f)
        print(mdf, file=f)
    f.close()

    import copy
    instruments = mdf.unstack().index
    iterables=[parameters_to_difference, instruments]
    new_columns = pd.MultiIndex.from_product(iterables, names=["var", "delta"])
    def converg_numofobs():
        d={}
        inst_list=[]
        for instrument, measurements in mdf.reset_index().groupby('pname'):
            inst_list.append(instrument)

        for instrument, measurements in mdf.reset_index().groupby('pname'):
            submdf=copy.deepcopy(measurements.xs(new_columns, axis=1))
            subcols= copy.deepcopy(submdf.columns)

            column_names_dict={}
            newsubd={}
            for colnum, colname in enumerate(subcols):
                column_name_entry=dict.fromkeys(['Param', 'Comparing_Inst'], {})
                column_name_entry['Param']=colname[0]
                column_name_entry['Comparing_Inst']=colname[1]
                column_names_dict[colnum]=column_name_entry

            for param in parameters_to_difference:
                instdict=dict.fromkeys(inst_list, {})
                for colnum in range(len(subcols)):
                    if column_names_dict[colnum]['Param'] == param:
                        subseries=submdf.iloc[:,colnum].dropna()
                        instdict[column_names_dict[colnum]['Comparing_Inst']]=subseries.expanding().mean().values
                newsubd[param]=instdict
                pprint.pprint(newsubd)

            d[instrument]=newsubd
        return d

    CON_NUMOB=converg_numofobs()
    #  pprint.pprint(CON_NUMOB)

    def singleday_stats(stat_type):
        d={}
        for instrument, measurements in mdf.reset_index().groupby('pname'):
            if stat_type == 'mean':
                d[instrument]  = measurements.xs(new_columns, axis=1).mean()
                table_lab_start='mean'
            elif stat_type == 'std':
                d[instrument]  = measurements.xs(new_columns, axis=1).std()
                table_lab_start='\n\nstandard deviation'
            elif stat_type == 'count':
                d[instrument]  = measurements.xs(new_columns, axis=1).count()
                table_lab_start='\n\nnumber of samples'
        return d, table_lab_start


    with open(outdir+'outputstats_'+day+'_'+timerange_str+'.txt', 'w') as f:
        for stat in stats_to_compute: #Multiday_df.index['stat']:# :
            d, lab_start=singleday_stats(stat)
            #  print(d.keys())
            print(lab_start+' of diffs between instrument measurements at approx the same time and loc', file=f)
            print(pd.DataFrame(data=d), file=f)
    f.close()
    # # # # # # # # # # # # # # # # # # # # # #
    for stat in stats_to_compute: #Multiday_df.index['stat']:# :
        d, _=singleday_stats(stat)
        for key in d.keys():
            TD= d[key].array
            T=d[key].index
            for entry in range(len(T)):
                Multiday_df.loc[(day,T[entry][0],stat), str(key)+'-'+str(T[entry][1])]=TD[entry]

    return Multiday_df, CON_NUMOB
# # # # # # # # # # # # # # # # # # # # # #
import matplotlib.ticker as ticker
def plot_CON_NUMOB(CON_NUMOB,day):
    #  fig, axs = plt.subplots(3)
    pprint.pprint(CON_NUMOB)
    print('\n')
    num_of_inst=len(CON_NUMOB.keys())
    #  for axnum, key in enumerate(CON_NUMOB):
    for first_inst in CON_NUMOB:
        fig, axs = plt.subplots(num_of_inst, figsize=(10,10))
        for axnum, second_inst in enumerate(CON_NUMOB[first_inst][parameters_to_difference[0]]):
            for param in parameters_to_difference:
                plotting_data = CON_NUMOB[first_inst][param][second_inst]

                if param == 'Thetav': C='blue'
                if param == 'Thetae': C='green'
                if param == 'tfast': C='red'


                axs[axnum].scatter(range(len(plotting_data)), plotting_data, color=C)
                if len(plotting_data) != 0:
                    axs[axnum].axhline(plotting_data[-1], linestyle=':',color=C, alpha=.5)
                axs[axnum].set_ylabel(second_inst)
                axs[axnum].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        axs[axnum].set_xlabel('count')

        plt.tight_layout()
        plt.savefig(plot_config.g_TORUS_directory+'/comparison/'+day+'_'+first_inst+'.png')
        plt.close()








# # # # # # # # # # # # # # # # # # # # # #

for day in day_list:
    #This is the directory path for the output file
    outdir = plot_config.g_TORUS_directory+day+'/data/intercomp/'
    # Setup Function Cache for speedup, if path directory does not make the directory
    if not os.path.exists(outdir): Path(outdir).mkdir(parents=True, exist_ok=True)

    Multiday_df, CON_NUMOB=approach_single_data_frame(day, plot_config, outdir, parameters_to_difference, Multiday_df)

    plot_CON_NUMOB(CON_NUMOB, day)

# * * * * *
for stat in stats_to_compute:
    for param in parameters_to_difference:
        Multiday_df.loc['Average', param, stat]=Multiday_df.loc[:,param,stat].mean(axis=0)
        Multiday_df.loc['Max', param, stat]=Multiday_df.loc[:,param,stat].max(axis=0)
        Multiday_df.loc['Min', param, stat]=Multiday_df.loc[:,param,stat].min(axis=0)



# # # # # # # # # # # # # # # # # # # # # #
with open(plot_config.g_TORUS_directory+'/comparison/stats/All_Multiday_stats.txt', 'w') as f:
    print(Multiday_df.T, file=f)
f.close()

for stat in stats_to_compute:
    with open(plot_config.g_TORUS_directory+'/comparison/stats/'+stat+'_Multiday_stats.txt', 'w') as f:
        if stat == 'count':
            print(Multiday_df.loc[:,parameters_to_difference[0],stat].T, file=f)
        else:
            print(Multiday_df.loc[:,:,stat].T, file=f)
    f.close()
for param in parameters_to_difference:
    with open(plot_config.g_TORUS_directory+'/comparison/stats/'+param+'_Multiday_stats.txt', 'w') as f:
        print(Multiday_df.loc[:,param,:].T, file=f)
    f.close()
for stat in stats_to_compute:
    if stat != 'count':
        for param in parameters_to_difference:
            with open(plot_config.g_TORUS_directory+'/comparison/stats/'+param+'_'+stat+'_Multiday_stats.txt', 'w') as f:
                print(Multiday_df.loc[:,param,stat].T, file=f)
                #  print(T_Multiday_df.loc[:,[:,param,:]], file=f)
            f.close()


for param in parameters_to_difference:
    meandf=Multiday_df.loc[:,param,'mean'].T
    stddf=Multiday_df.loc[:,param,'std'].T
    countdf=Multiday_df.loc[:,param,'count'].T
    for index,meanrow in meandf.iterrows():

        icombo=index
        means=meanrow[:-3].to_numpy()
        std_interm=stddf.loc[index]
        std=std_interm[:-3].to_numpy()
        count_interm=countdf.loc[index]
        samples=count_interm[:-3].to_numpy()
        
        fig, axs = plt.subplots(2)
        axs[0].scatter(samples, means)
        axs[0].set_xlabel('count')
        axs[0].set_ylabel('bias mean')
        axs[1].scatter(samples, std, color='red')
        axs[1].set_xlabel('count')
        axs[1].set_ylabel('bias std')
        plt.title(icombo)
        plt.tight_layout()
        plt.savefig(plot_config.g_TORUS_directory+'/comparison/stats/images/'+param+'_'+icombo+'.png')
        plt.close()


exit()
