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


height_level = 17 #0.81 eta level
height_level = 3 #roughly 80 m above ground level
grid_spacing = 12*1000 #km 2 m

xdim = 405
ydim = 325
time_step = 21 #2100hrs

figwidth = 6
FigSize=(figwidth, ydim/xdim*figwidth)

for species in ['SO2','ANAJ','NO2','O3']:
    tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
    ncfile="species_data_"+species+".nc"#"ftle_80m.nc"
    root = Dataset(ncfile,'r')
    vars = root.variables
    u = vars['eastward_vel'][:]
    v = vars['northward_vel'][:]
    lat = vars['lat'][:]
    lon = vars['lon'][:]
    root.close()
    u=u[time_step,:,:]
    v=v[time_step,:,:]
    
    
    #lon = lon[-1,:,:]
    #lat = lat[-1,:,:]
    #latin = lat[:25,:,:]
    #longin = lon[:25,:,:]
    
    #x = np.linspace(-grid_spacing*(xdim-1)/2,grid_spacing*(xdim-1)/2,xdim)
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
    cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True)
    #cs=m.contourf(lon,lat,s1,levels=np.linspace(s1.min(axis=None),s1.max(axis=None),301),latlon=True,vmin=0.5*s1.min(axis=None),vmax=0.5*s1.max(axis=None))
    m.drawcoastlines()
    m.drawstates()
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-2)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    t=1
    hrs, mins = np.divmod((t-1)*10,60)
    plt.title("Integration time = -{0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
    plt.savefig(species+'_lcs_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
    plt.close('all')
