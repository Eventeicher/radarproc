import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

#  %pylab inline
import glob
#  import mpl_toolkits
#  print(mpl_toolkits)
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm


def sparseFull(a):

    """This expands a sparse grid, incorporating all data values"""

    # Create empty shape for filling
    b = np.full([len(a.dimensions['Lat']),len(a.dimensions['Lon'])], -99900.)

    # Create empty shape for filling
    xlen = len(a.dimensions['Lon'])
    ylen = len(a.dimensions['Lat'])

    # Gather into four arrays
    varname = a.TypeName
    var = a.variables[varname][:]
    pixel_x = a.variables['pixel_x'][:]
    pixel_y = a.variables['pixel_y'][:]
    pixel_count = a.variables['pixel_count'][:]

    # Loop through the four arrays simultaneously, populating 2D array 'b'
    for w,x,y,z in zip(var, pixel_x, pixel_y, pixel_count):
        distToEdge = xlen - y
        # check if pixel count exceeds distance to right edge of domain
        if z > distToEdge:
            b[x,y:xlen] = w
            pixelsLeft = z - distToEdge
            # check if pixels remaining are less than grid width
            if pixelsLeft <= xlen:
                b[x+1,0:(pixelsLeft-1)] = w
            else:
                rowsLeft, pixelCount_RemainingRow = divmod(pixelsLeft, xlen)
                b[x:(x+rowsLeft+1),0:7001] = w
                # check if pixels are remaining
                if pixelCount_RemainingRow > 0:
                    b[(x+rowsLeft+1),0:pixelCount_RemainingRow] = w
        else:
            b[x,y:(y + z)] = w

    return b

def sparseReal(a, acc_method):
    """This expands a sparse grid over a -99900 or "0" grid, ignoring all -999xx values
       from the sparse grid
    FAST"""
    # Size of grid from netCDF file
    xlen = a.dimensions['Lon'].size
    ylen = a.dimensions['Lat'].size
    # Create empty shape for filling
    if acc_method == 1:
        b = np.full([ylen, xlen], -99900.)
    else:
        b = np.zeros([ylen, xlen])
    return a, b

def readwdssii(fin):

    """This loads in the file"""

    a = nc.Dataset(fin, 'r')

    xlen = len(a.dimensions['Lon'])
    ylen = len(a.dimensions['Lat'])
    lat = a.Latitude
    lon = a.Longitude
    try:
        varname = a.TypeName


        var = a.variables[varname][:]

        if a.DataType[0:6] == 'Sparse':

            #var = sparseFull(a)
            var = sparseReal(var,b)

    except:
        print('except')
        return xlen,ylen,varname,0
    a.close()
    del(a)

    return xlen, lat, ylen, lon, varname, var


###CODE HERE
####can loop over this...
path = '/home/chris/codes/data/20190608/merged/MergedReflectivityQC/00.25/20190609-010629.netcdf'
#  path = '/home/chris/codes/data/20190608/merged/MergedAzShear_0-2kmAGL/00.00/20190609-010746.netcdf'

gridspacing=.01
#  gridspacing=.005
x1,ulat,y1,ulon, varname,ref1 = readwdssii(path)
#  x1,ulat,y1,ulon, varname,ref1 = readwdssii('/Users/aereinha/Downloads/RotationTracks60min/00.50/20170528-235400.netcdf')


#longitude runs west to east
#latitude runs north to south
lat = np.linspace(ulat,stop=ulat+(y1*-gridspacing),num=y1)
lon = np.linspace(ulon,stop=ulon+(x1*gridspacing),num=x1)
#  lat = np.linspace(ulat,stop=ulat+(y1*-0.01),num=y1)
#  lon = np.linspace(ulon,stop=ulon+(x1*0.01),num=x1)
nx,ny = np.meshgrid(lon,lat)
fig = plt.figure()
#ulat is upper corner, lon[-1] is lower corner lat[-1] is upper and ulat is lower
m = Basemap(resolution = 'l',llcrnrlat=lat[-1],llcrnrlon=ulon,urcrnrlat=ulat,urcrnrlon=lon[-1],projection='lcc',lat_0=38.0,lon_0=-97)

xx,yy = m(nx,ny)
#need to divide AzShear by 1000 because of scaling...
m.contourf(xx,yy,ref1,levels=np.arange(0,70,10),cmap=plt.get_cmap('viridis'))
plt.colorbar()
m.drawstates()
m.drawcountries()
#  m.drawcoastlines()
#  plt.savefig()
plt.show()
