#!/usr/bin/env python

# -*- coding: utf-8 -*-
from joblib import Memory
import glob
import pickle
import os
import fire
import pyart

import pathlib

###Set Up Temporary Directory
############################
cachedir='./cachedir'
mem= Memory(cachedir,verbose=1)
how_many_fail= 0

def dealias_radar(radar_pickled):

    # Regenerate the calculated values in the object
    radar = pickle.loads(radar_pickled)
    dealias_data = pyart.correct.region_dealias.dealias_region_based(radar)

    return dealias_data


def info_radar_files(level, files):
    for radar_file in files[:]:

        print("\nInfo ... " +
              " Input File: " + str(radar_file))

        radar = pyart.io.read(radar_file)

        radar.info(level = level)


def dealias_radar_files(prefix, files):

    #joblib dealias stuff
    persistant_dealias_radar = mem.cache( dealias_radar )

    for radar_file in files[:]:
        #print("a " + radar_file)
        #print("b " + str(radar_file))
        #print("c " + os.path.dirname(radar_file))
        #print("d " + prefix)
        #print("e " + os.path.basename(radar_file))

        out_filename = os.path.dirname(radar_file) + '/' + prefix + os.path.basename(radar_file)

        print("\nDealiasing... " +
              " Input File: " + str(radar_file) +
              " Output File: " + out_filename)

        if radar_file == '.':
            how_many_fail = how_many_fail + 1
            continue;

        radar = pyart.io.read(radar_file)

        if radar.scan_type == 'rhi':
            print('Skipping, we are not dealiasing RHI files at this time \n')
            continue;

        try:
            if radar.iter_field('corrected_velocity'):
                print('Input file already contains dealias data in field corrected_velocity...... so just creating a copy with the specified prefix \n')
                pyart.io.write_cfradial(out_filename, radar)
                continue;
        except:
            print('')

        dealias_data = persistant_dealias_radar( pickle.dumps(radar) )
        radar.add_field('corrected_velocity', dealias_data, replace_existing=True)

        pyart.io.write_cfradial(out_filename, radar)

class Radar_File_Functions(object):

    def dealias_files(self, *input_files, output_file_prefix: str = 'dealiased_'):
        '''
        Dealias radar file(s).

        input_files : array
          Radar files to process.
        output_file_prefix : string
          Path prefix to start all output files with.
        '''
        dealias_radar_files(output_file_prefix, input_files)

    def dump_file_info(self, *input_files, pyart_info_level: str = 'standard'):
        '''
        Describe radar file(s).

        input_files : array
          Radar files to process.
        pyart_info_level : string
          Pyart specification for info level.
        '''
        info_radar_files(pyart_info_level, input_files)

def main():
  fire.Fire(Radar_File_Functions)

if __name__ == '__main__':
    main()

print('*****************\n', how_many_fail,'\n***********************')
