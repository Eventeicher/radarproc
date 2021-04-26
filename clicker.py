#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import needed modules
######################
import matplotlib
#  matplotlib.use('agg')
#  matplotlib.use('GTK3Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.patheffects as PE
from matplotlib.collections import PatchCollection
from matplotlib.offsetbox import (AnchoredOffsetbox, DrawingArea, HPacker, TextArea)
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredAuxTransformBox
from matplotlib import ticker
from matplotlib.ticker import (LinearLocator, FixedLocator, MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.colors 
import time
from cycler import cycler
from mpl_toolkits.axes_grid1.axes_divider import make_axes_area_auto_adjustable
import matplotlib.cm as cmx
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import LineCollection
from cycler import cycler
from datetime import datetime, timedelta
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
import pandas as pd
from pandas.plotting import register_matplotlib_converters
import numpy as np
import xarray as xr
from netCDF4 import num2date
import logging
import csv 
import os
import os.path
import time
from pathlib import Path
from joblib import Memory, Parallel, delayed
from scipy import ndimage, interpolate
import cmocean
import gc
import glob
import nexradaws
import pyart
import gc
from pympler.asizeof import asizeof
#this is the file with the plotting controls to access any of the vars in that file use config.var
import config as plot_config
import copy
import singledop
import pynimbus as pyn 
import tracemalloc
import shapefile
import math 
#  tracemalloc.start()

import pprint

def clicker_defn(plt):
    pts = [np.nan]
    #  plt.setp(plt.gca(), autoscale_on=False)
    #  plt.draw()
    #  plt.waitforbuttonpress()
    pts = plt.ginput(n=40,timeout=0) # look up ginput docs for some good guidance
    #  plt.waitforbuttonpress()
    return pts

def make_csv(Data, pts, surge_name, filename, Does_csv_exist):
        if Does_csv_exist == False: 
            # making csv file if it does not exist 
            with open(filename, 'w') as csvfile: 
                # field names 
                fields = ['RName', 'Scan_time','Surge_ID', 'points'] 
                
                # creating a csv writer object 
                csvwriter = csv.writer(csvfile) 
                    
                # writing the fields 
                csvwriter.writerow(fields) 
        ###
        # Now append the new row of data points
        with open(filename, 'a') as csvfile: 
            # writing the data rows 
            csvwriter = csv.writer(csvfile) 
            rows= [[Data['P_Radar'].name,Data['P_Radar'].Scan_time, surge_name, pts]]
            csvwriter.writerows(rows)
