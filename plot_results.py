from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
from scipy import interpolate
import mapping_functions as mf
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker as ticker
import time
import calendar
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}
plt.close('all')

xdim = 405
ydim = 325
tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
#ncfile="vapor_flux_ftle.nc"#"ftle_80m.nc"
species = 'NA'
ncfile=species+"_FTLE.nc"#"ftle_80m.nc"
#ncfile="ftle_vapor_flux.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
'''
lat = vars["initial_lat"][:].reshape([ydim,xdim])
lon = vars["initial_lon"][:].reshape([ydim,xdim])
'''
lat = vars["lat"][:]#.reshape([ydim,xdim])
lon = vars["lon"][:]#.reshape([ydim,xdim])
print(vars['FTLE'][:].shape)
time = 86400*vars["time"][:]+tstart
tdim = np.shape(time)[0]
ftle = vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
#lon, lat = np.meshgrid(lon,lat)
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
            resolution = 'h',
            area_thresh=1000.,
            )


lon,lat = np.meshgrid(lon,lat)
cs=m.contourf(lon,lat,ftle[-1,:,:],levels=np.linspace(ftle.min(axis=None),ftle.max(axis=None),301),latlon=True)
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#cs=m.pcolor(lon,lat,ftle[-1,:,:],latlon=True)#,levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),301),latlon=True)
for t in range(2,tdim+1):#146):#time)):
    for c in cs.collections:
        c.remove()
    #cs.set_array(np.ravel(ftle[-t,:,:]))
    #cs=m.contourf(lon,lat,ftle[:,:,-t],levels=np.linspace(ftle.min(axis=None),ftle.max(axis=None),301),latlon=True)
    print(t)
    cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),301),latlon=True,vmin=0.5*ftle[-t,:,:].min(axis=None),vmax=0.5*ftle[-t,:,:].max(axis=None))
    #cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),301),latlon=True)
    hrs, mins = np.divmod((t-1)*10,60)
    plt.title("Integration time = -{0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
    plt.savefig(species+'_lcs_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
#'''