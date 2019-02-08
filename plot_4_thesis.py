from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
plt.close('all')
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}


height_level = 17 #0.81 eta level
height_level = 3 #roughly 80 m above ground level
grid_spacing = 12*1000 #km 2 m

xdim = 405
ydim = 325
time_step = 21 #2100hrs

figwidth = 6.5
FigSize=(figwidth, ydim/xdim*figwidth)
#FigSize=(3*figwidth, 8*figwidth)
fig = plt.figure(1,figsize=FigSize)
ncfile="species_data_SO2.nc"#"ftle_80m.nc"
root = Dataset(ncfile,'r')
vars = root.variables
lat = vars['lat'][:]
lon = vars['lon'][:]
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
            suppress_ticks=True,
            )



ftle_species = ['wind_speed','water_vapor','NO2','O3']


ncfile=ftle_species[0]+"_FTLE.nc"#"ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
lat = vars["lat"][:]
lon = vars["lon"][:]
ftle = vars['FTLE'][:,:,:,0]
lon, lat = np.meshgrid(lon,lat)
root.close()
plt.subplot(221)
cs=m.contourf(lon,lat,ftle[-7,:,:],levels=np.linspace(np.min(ftle[-7,:,:],axis=None),np.max(ftle[-7,:,:],axis=None),301),latlon=True)
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-3)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
plt.title("Horizontal Wind",**titlefont)


ncfile=ftle_species[1]+"_FTLE.nc"#"ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
lat = vars["lat"][:]
lon = vars["lon"][:]
ftle = vars['FTLE'][:,:,:,0]
lon, lat = np.meshgrid(lon,lat)
root.close()
plt.subplot(222)
cs=m.contourf(lon,lat,ftle[-7,:,:],levels=np.linspace(np.min(ftle[-7,:,:],axis=None),np.max(ftle[-7,:,:],axis=None),301),latlon=True)
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-3)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
plt.title("Water Vapor",**titlefont)

ncfile=ftle_species[2]+"_FTLE.nc"#"ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
lat = vars["lat"][:]
lon = vars["lon"][:]
ftle = vars['FTLE'][:,:,:,0]
lon, lat = np.meshgrid(lon,lat)
root.close()
plt.subplot(223)
cs=m.contourf(lon,lat,ftle[-7,:,:],levels=np.linspace(np.min(ftle[-7,:,:],axis=None),np.max(ftle[-7,:,:],axis=None),301),latlon=True)
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-3)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
plt.title("Nitrogen Dioxide",**titlefont)

ncfile=ftle_species[3]+"_FTLE.nc"#"ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
lat = vars["lat"][:]
lon = vars["lon"][:]
ftle = vars['FTLE'][:,:,:,0]
lon, lat = np.meshgrid(lon,lat)
root.close()
plt.subplot(224)
cs=m.contourf(lon,lat,ftle[-7,:,:],levels=np.linspace(np.min(ftle[-7,:,:],axis=None),np.max(ftle[-7,:,:],axis=None),301),latlon=True)
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-3)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
plt.title("Ozone",**titlefont)

plt.savefig('pollution.png', transparent=False, bbox_inches='tight',dpi=200)





