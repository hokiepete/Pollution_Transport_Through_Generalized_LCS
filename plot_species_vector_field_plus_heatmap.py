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
#,'ANH4J'
xdim = 405
ydim = 325
figwidth = 6.5
#FigSize=(figwidth, 2*ydim/xdim*figwidth)
#subs = [421,422,423,424,425,426,427]
FigSize=(figwidth, 0.5*ydim/xdim*figwidth)
subs = [241,242,243,244,245,246,247]
for s in ['NO2','SO2','O3','ASO4J','ANAJ']:
    fig = plt.figure(1,figsize=FigSize)
    for k, time in enumerate([21,20,19,18,17,16,15]):
        ax = plt.subplot(subs[k])
        root = nc.Dataset('species_data_'+s+'.nc','r')
        vars = root.variables
        species = vars['species'][time,:,:]
        u = vars['eastward_vel'][time,:,:]
        v = vars['northward_vel'][time,:,:]
        lon = vars['lon'][:]
        lat = vars['lat'][:]
        
        lon_min = np.min(lon,axis=None)
        lon_max = np.max(lon,axis=None)
        lat_min = np.min(lat,axis=None)
        lat_max = np.max(lat,axis=None)
        parallels = np.arange(round(lat_min,0),lat_max+2,2)
        meridians = np.arange(round(lon_max,0),lon_min-2,-3)
        
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
        #m.contourf(lon,lat,species,levels=np.linspace(np.min(species,axis=None),np.max(species,axis=None),301),latlon=True,cmap=matplotlib.cm.jet)
        m.pcolormesh(lon,lat,species,latlon=True,cmap=matplotlib.cm.jet,shading='gouraud')
        m.colorbar(format='%.3f')
        #cb.set_ticklabels(fmt)
        #if k == 0:
        ax.set_title("Time = {:d}:00 hrs".format(time))
        #else:
        #    ax.set_title("Time = {:d}".format(-k))
        velquiver = m.quiver(lon[::downsamp,::downsamp],lat[::downsamp,::downsamp],u[::downsamp,::downsamp],v[::downsamp,::downsamp],latlon=True)
        qk = plt.quiverkey(velquiver, 0.9*m.urcrnrx, 1.03*m.urcrnry, 10, '$10 $', labelpos='E', coordinates='data', fontproperties={'size': '10'})
        m.drawcoastlines()
        m.drawstates()
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
    #plt.tight_layout()
    plt.subplots_adjust(wspace=1.0,hspace=0.3)
    #plt.savefig(s+'quiver_vert.png',dpi=300)
    plt.savefig(s+'quiver_hori.png',dpi=300)
    #plt.savefig(s+'quiver.eps')
    plt.close('all')
























