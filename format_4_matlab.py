from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
#from scipy import interpolate
import mapping_functions as mf
#from scipy.io import savemat
from hdf5storage import savemat
height_level = 17 #0.81 eta level
height_level = 3 #roughly 80 m above ground level
grid_spacing = 12 #km

tdim = 25
xdim = 102
ydim = 82

root = Dataset('wrf_2011_07_01','r')
cen_lat = getattr(root,'CEN_LAT')
cen_lon = getattr(root,'CEN_LON')
true_lat1 = getattr(root,'TRUELAT1')
true_lat2 = getattr(root,'TRUELAT2')
ref_lat = getattr(root,'MOAD_CEN_LAT')
ref_lon = getattr(root,'STAND_LON')
vars = root.variables
#Wind Velocity
u = mf.unstagger(vars['U'][:],axis=3)
v = mf.unstagger(vars['V'][:],axis=2)
w = mf.unstagger(vars['W'][:],axis=1)
grd = vars['HGT'][0,:,:]
hgt = mf.unstagger((vars['PH'][0,:,:,:] + vars['PHB'][0,:,:,:])/9.81,axis=0)-grd
#Water Vapor Flux, Vertically Integrated
lat_in = vars['XLAT'][0,:,:]
lon_in = vars['XLONG'][0,:,:]
#root.close()

xin = np.linspace(0,12000*101,102)
yin = np.linspace(0,12000*81,82)
hin = hgt.mean(axis=(1,2))
time = np.linspace(0,3600*23,24)
u = u.data
v = v.data
w = w.data
hin = hin.data

savemat('wrfdata.mat',{'x':xin,'y':yin,'z':hin,'time':time,'u':u,'v':v,'w':w})


"""
xin,yin = mf.lonlat2m(ref_lon,ref_lat,lon_in,lat_in,true_lat1,true_lat2)
xin=xin-xin.min()
yin=yin-yin.min()
xin = np.linspace(0,12000*101,102)
x = np.empty([35,ydim,xdim])
y = np.empty([35,ydim,xdim])
for h in range(35):
    x[h,:,:]=xin
    y[h,:,:]=yin

time = np.linspace(0,3600*23,24)
"""