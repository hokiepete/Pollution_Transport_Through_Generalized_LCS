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

time = 0
root = nc.Dataset('redux_species_data_Q'+'.nc','r')
vars = root.variables
lon = vars['lon'][:]
lat = vars['lat'][:]
root.close()
lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)

parallels = [31,34,37]#np.arange(round(lat_min,0),lat_max+2,2)
meridians = [-86,-82,-78]#np.arange(round(lon_max,0),lon_min-2,-3)

m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'l',
            area_thresh=1000.,
            suppress_ticks=True,
            )

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}
time = 21
s = 'O3'
#,'ANH4J'
xdim = 405
ydim = 325
figwidth = 7
#FigSize=(figwidth, 2*ydim/xdim*figwidth)
#subs = [421,422,423,424,425,426,427]
FigSize=(figwidth, .525*ydim/xdim*figwidth)
subs = [231,232,233,234,235,236]
label = ['H2O', 'SO2', 'O3', 'SO4', 'NO2', 'Na']
for time in [15]:
    fig = plt.figure(1,figsize=FigSize)
    for k, s in enumerate(['Q','SO2','O3','ASO4J','NO2','ANAJ']):
        print(k,s,subs[k])
        ax = plt.subplot(subs[k])
        root = nc.Dataset('redux_species_data_'+s+'.nc','r')
        vars = root.variables
        if s=='Q':
            species = vars['species'][time,:,:]
            u = vars['eastward_vel'][time,:,:]
            v = vars['northward_vel'][time,:,:]
        
        else:
            species = vars['species'][time,:,:]/1000
            u = vars['eastward_vel'][time,:,:]
            v = vars['northward_vel'][time,:,:]
        
        lon = vars['lon'][:]
        lat = vars['lat'][:]
        
        """
        lon_min = np.min(lon,axis=None)
        lon_max = np.max(lon,axis=None)
        lat_min = np.min(lat,axis=None)
        lat_max = np.max(lat,axis=None)
        parallels = [31,34,37]#np.arange(round(lat_min,0),lat_max+2,2)
        meridians = [-86,-82,-78]#np.arange(round(lon_max,0),lon_min-2,-3)
        
        m = Basemap(llcrnrlon=lon_min,
                    llcrnrlat=lat_min,
                    urcrnrlon=lon_max,
                    urcrnrlat=lat_max,
                    projection='merc',
                    resolution = 'l',
                    area_thresh=1000.,
                    suppress_ticks=True,
                    )
        """
        downsamp = 25
        lon,lat = np.meshgrid(lon,lat)
        print(lon.shape)
        #m.contourf(lon,lat,species,levels=np.linspace(np.min(species,axis=None),np.max(species,axis=None),301),latlon=True,cmap=matplotlib.cm.jet)
        #m.pcolormesh(lon,lat,species,latlon=True,cmap=matplotlib.cm.jet,shading='gouraud')
        m.pcolormesh(lon,lat,species,latlon=True,cmap='Greys',shading='gouraud')
        plt.annotate('{0}'.format(label[k]),xy=(0.03*m.urcrnrx,0.03*m.urcrnry),fontsize=7)
        #
        #cb.set_ticklabels(fmt)
        if k == 4:
            cbar = m.colorbar(format='%.3f')
        else:
            cbar = m.colorbar()
            
        cbar.ax.tick_params(labelsize=6) 
        tick_locator = ticker.MaxNLocator(nbins=9)
        cbar.locator = tick_locator
        cbar.update_ticks()
        
        #ax.set_title("Time = {:d}:00 hrs".format(time))
        #else:
        #    ax.set_title("Time = {:d}".format(-k))
        #velquiver = m.quiver(lon[::downsamp,::downsamp],lat[::downsamp,::downsamp],u[::downsamp,::downsamp],v[::downsamp,::downsamp],latlon=True,color='gray')
        velquiver = m.quiver(lon[::downsamp,::downsamp],lat[::downsamp,::downsamp],u[::downsamp,::downsamp],v[::downsamp,::downsamp],latlon=True,color='k')
        #if k in [1,4]:
        qk = plt.quiverkey(velquiver, 0.75*m.urcrnrx, 1.05*m.urcrnry, 10, '$10 m/s$', labelpos='E', coordinates='data', fontproperties={'size': '6'})
        #else:
        #    qk = plt.quiverkey(velquiver, 0.9*m.urcrnrx, 1.03*m.urcrnry, 100, '$100 m/s$', labelpos='E', coordinates='data', fontproperties={'size': '6'})

        m.drawcoastlines()
        m.drawstates()
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=7)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=7)
    #plt.tight_layout()
    #plt.subplots_adjust(wspace=1.0,hspace=0.3)
    #plt.savefig(s+'quiver_vert.png',dpi=300)
    plt.savefig('quiver_comparison.png',dpi=300)
    #plt.savefig(s+'quiver.eps')
    #plt.close('all')
























