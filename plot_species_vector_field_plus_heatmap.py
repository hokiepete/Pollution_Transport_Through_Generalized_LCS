import netCDF4 as nc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
os.environ["PROJ_LIB"] = "C:\\ProgramData\\Anaconda3\\Library\\share"
from mpl_toolkits.basemap import Basemap

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \cdot 10^{{{}}}$'.format(a, b)


plt.close('all')
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}
time = 21
s = 'O3'
#for s in ['NO2','SO2','O3','ANH4J','ASO4J','ANAJ']:
for time in [21,20,19,18,17,16,15]:
    root = nc.Dataset('species_data_'+s+'.nc','r')
    vars = root.variables
    species = vars['species'][time,:,:]
    u = vars['eastward_vel'][time,:,:]
    v = vars['northward_vel'][time,:,:]
    lon = vars['lon'][:]
    lat = vars['lat'][:]
    xdim = 405
    ydim = 325
    
    lon_min = np.min(lon,axis=None)
    lon_max = np.max(lon,axis=None)
    lat_min = np.min(lat,axis=None)
    lat_max = np.max(lat,axis=None)
    
    figwidth = 6.5
    FigSize=(figwidth, ydim/xdim*figwidth+0.5)
    #FigSize=(3*figwidth, 8*figwidth)
    fig = plt.figure(1,figsize=FigSize)
    
    m = Basemap(llcrnrlon=lon_min,
                llcrnrlat=lat_min,
                urcrnrlon=lon_max,
                urcrnrlat=lat_max,
                projection='merc',
                resolution = 'l',
                area_thresh=1000.,
                suppress_ticks=True,
                )
    
    downsamp = 25
    lon,lat = np.meshgrid(lon,lat)
    m.contourf(lon,lat,species,levels=np.linspace(np.min(species,axis=None),np.max(species,axis=None),301),latlon=True)
    #m.colorbar(format=ticker.FuncFormatter(fmt))
    m.colorbar(format='%.1e')
    #cb.set_ticklabels(fmt)
    velquiver = m.quiver(lon[::downsamp,::downsamp],lat[::downsamp,::downsamp],u[::downsamp,::downsamp],v[::downsamp,::downsamp],latlon=True)
    qk = plt.quiverkey(velquiver, 0.9*m.urcrnrx, 1.03*m.urcrnry, 10, '$10 $', labelpos='E', coordinates='data', fontproperties={'size': '10'})
    m.drawcoastlines()
    m.drawstates()
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-3)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
    plt.savefig(s+'quiver_{:d}.png'.format(time),padding=3)
    plt.close('all')
























