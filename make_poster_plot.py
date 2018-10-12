from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
from scipy import interpolate
import mapping_functions as mf
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import time
import calendar
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

figwidth = 4
FigSize=(2.75*figwidth, 8*ydim/xdim*figwidth)
#FigSize=(3*figwidth, 8*figwidth)
fig = plt.figure(1,figsize=FigSize)
gs1 = gridspec.GridSpec(8, 3)
gs1.update(wspace=0.05, hspace=0.05)
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



ftle_species = ['wind_speed','water_vapor','SO2','NA','NO2','O3','SO4','NH4']
final_len = len(ftle_species)-1
for k, species in enumerate(['wind_speed','water_vapor','SO2','ANAJ','NO2','O3','ASO4J','ANH4J']):
    tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
    ncfile="species_data_"+species+".nc"#"ftle_80m.nc"
    root = Dataset(ncfile,'r')
    vars = root.variables
    u = vars['eastward_vel'][:]
    v = vars['northward_vel'][:]
    lat = vars['lat'][:]
    lon = vars['lon'][:]
    lon,lat = np.meshgrid(lon,lat)
    root.close()
    u=u[time_step,:,:]
    v=v[time_step,:,:]
    
    
    x = np.linspace(0,grid_spacing*(xdim-1),xdim)
    dx = x[1]-x[0]
    #y = np.linspace(-grid_spacing*(ydim-1)/2,grid_spacing*(ydim-1)/2,ydim)
    y = np.linspace(0,grid_spacing*(ydim-1),ydim)
    dy = y[1]-y[0]
    x, y = np.meshgrid(x,y)
    
    dudy,dudx = np.gradient(u,dy,dx)
    dvdy,dvdx = np.gradient(v,dy,dx)
    
    
    s1 = np.ma.empty([ydim,xdim])
    J = np.array([[0, 1], [-1, 0]])
    for i in range(ydim):
        for j in range(xdim):
            if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
                Utemp = np.array([u[i, j], v[i, j]])
                Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
                S = 0.5*(Grad + np.transpose(Grad))
                s1[i,j] = -3600*np.min(np.linalg.eig(S)[0])
    
            else:
                s1[i,j] = np.ma.masked
        del j
    del i
    
    #plt.subplot(8,3,int(3*k+1))
    plt.subplot(gs1[k,0])
    cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True)
    #cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
    m.drawcoastlines()
    m.drawstates()
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-2)
    #m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    #m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    plt.xticks([])
    plt.yticks([])
    
    ncfile=ftle_species[k]+"_FTLE.nc"#"ftle_80m.nc"
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    lat = vars["lat"][:]
    lon = vars["lon"][:]
    ftle = vars['FTLE'][:,:,:,0]
    lon, lat = np.meshgrid(lon,lat)
    root.close()

    #plt.subplot(8,3,int(3*k+2))
    plt.subplot(gs1[k,1])
    cs=m.contourf(lon,lat,ftle[-7,:,:],levels=np.linspace(np.min(ftle[-7,:,:],axis=None),np.max(ftle[-7,:,:],axis=None),301),latlon=True)
    #cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
    m.drawcoastlines()
    m.drawstates()
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-2)
    #m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    #m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    plt.xticks([])
    plt.yticks([])
    
    #plt.subplot(8,3,int(3*k+3))
    plt.subplot(gs1[k,2])
    cs=m.contourf(lon,lat,ftle[-37,:,:],levels=np.linspace(np.min(ftle[-37,:,:],axis=None),np.max(ftle[-37,:,:],axis=None),301),latlon=True)
    #cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
    m.drawcoastlines()
    m.drawstates()
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-2)
    #m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    #m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    plt.xticks([])
    plt.yticks([])

    print(k/final_len*100)

plt.savefig('poster_plot.tif', transparent=False, bbox_inches='tight',dpi=200)





