
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import calendar
import time
from mpl_toolkits.basemap import Basemap
plt.close('all')

xdim = 405
ydim = 325
tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
ncfile="ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
lat = vars["lat"][:]#.reshape([ydim,xdim])
lon = vars["lon"][:]#.reshape([ydim,xdim])
print(vars['FTLE'][:].shape)
time = 86400*vars["time"][:]+tstart
tdim = np.shape(time)[0]
#ftle = vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
root.close()
lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)

m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'c',
            area_thresh=1000.,
            )

#lon,lat = np.meshgrid(lon,lat)
#cs=m.contourf(lon,lat,ftle[-1,:,:],levels=np.linspace(ftle.min(axis=None),ftle.max(axis=None),301),latlon=True)
#ncfile="SE_tracers.nc"
ncfile="tracers.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
t_lat = vars["lat"][:]#.reshape([ydim,xdim])
t_lon = vars["lon"][:]#.reshape([ydim,xdim])
time = 86400*vars["time"][:]+tstart
tdim = np.shape(time)[0]
root.close()
attracting_ridges_lat=[]
attracting_ridges_lon=[]
for a in [13,17,40,46,47,55]:
    ncfile="attracting_ridge_{:02d}.nc".format(a)
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    attracting_ridges_lat.append(vars["lat"][:])
    attracting_ridges_lon.append(vars["lon"][:])
    root.close()

repelling_ridges_lat=[]
repelling_ridges_lon=[]
for r in [8,18,19,53,56,60]:
    ncfile="repelling_ridge_{:02d}.nc".format(r)
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    repelling_ridges_lat.append(vars["lat"][:])
    repelling_ridges_lon.append(vars["lon"][:])
    root.close()
for t in range(146):#time)):
    #for c in cs.collections:
    #    c.remove()
    #cs.set_array(np.ravel(ftle[:,:,t]))
    #cs=m.contourf(lon,lat,ftle[:,:,-t],levels=np.linspace(ftle.min(axis=None),ftle.max(axis=None),301),latlon=True)

    #cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),301),latlon=True)
    #plt.title("{0}".format(time[-t]),fontsize=18)
    plt.figure(figsize=[16,12])
    for i in range(len(attracting_ridges_lat)):
        m.plot(attracting_ridges_lon[i][:,t],attracting_ridges_lat[i][:,t],'b-',latlon=True)
    for i in range(len(repelling_ridges_lat)):
        m.plot(repelling_ridges_lon[i][:,t],repelling_ridges_lat[i][:,t],'r-',latlon=True)
    m.scatter(t_lon[:,t],t_lat[:,t],latlon=True)
    
    
    m.drawcoastlines()
    m.drawstates()
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-2)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)    
    plt.savefig('OECS_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
    plt.close('all')
#'''
    
    
    
"""
ncfile="SE_ridge.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
r_lat = vars["lat"][:]#.reshape([ydim,xdim])
r_lon = vars["lon"][:]#.reshape([ydim,xdim])
root.close()

"""