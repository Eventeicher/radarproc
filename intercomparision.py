#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.patheffects as PE
from matplotlib.collections import PatchCollection
from matplotlib import ticker
from matplotlib.ticker import (LinearLocator, FixedLocator, MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import ListedColormap, Normalize
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
import datetime as dt
from datetime import datetime, date, timedelta
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
import argparse, cProfile, logging, time, os, os.path
from os.path import expanduser
from pathlib import Path
from joblib import Memory, Parallel, delayed
from scipy import ndimage, interpolate
from operator import attrgetter
from collections import namedtuple
import pyart, nexradaws, sys, traceback, shutil, glob, gc, cmocean
import fnmatch
import gc
from pympler.asizeof import asizeof
#this is the file with the plotting controls to access any of the vars in that file use config.var
import config as plot_config
import copy

from read_pforms import read_TInsitu, pform_names

pd.set_option('display.max_columns', None)
df_holder =[]
p_holder=[]
for p in pform_names('TInsitu'):
    print(p)
    print('............')
    d_avail=read_TInsitu(plot_config, p, True)
    if d_avail == True:
        df, ptype =read_TInsitu(plot_config, p)
        df['pname'] = p
        df = df.drop(columns=['U','V', 'dir', 'spd' ])
        if ptype =='NSSL':
            df = df.drop(columns=['id', 'qc1', 'qc2', 'qc3', 'qc4'])
        df_holder.append(df)
        p_holder.append(p)
        #  print(df)
full_df= pd.concat(df_holder)
new_cols=['dist_Prb1', 'dist_Prb2', 'dist_FFld', 'dist_Wins', 'dist_LiDR', 'dist_CoMeT1', 'dist_CoMeT2']
for c in new_cols:
    full_df.insert(2, c, np.nan)
full_df.set_index(['pname', 'datetime'], inplace=True)
print(full_df)


