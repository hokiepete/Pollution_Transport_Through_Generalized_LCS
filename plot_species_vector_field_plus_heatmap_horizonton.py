import netCDF4 as nc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
os.environ["PROJ_LIB"] = "C:\\ProgramData\\Anaconda3\\Library\\share"
from mpl_toolkits.basemap import Basemap


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
    
    root = nc.Dataset('species_data_'+s+'.nc','r')
    vars = root.variables
    species_in = vars['species'][:]
    u_in = vars['eastward_vel'][:]
    v_in = vars['northward_vel'][:]
    lon = vars['lon'][:]
    lat = vars['lat'][:]
    species_max = species_in[15:22,:,:].max()
    species_min = species_in[15:22,:,:].min()
    lon_min = np.min(lon,axis=None)
    lon_max = np.max(lon,axis=None)
    lat_min = np.min(lat,axis=None)
    lat_max = np.max(lat,axis=None)
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-4)
    
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
    
    for k, time in enumerate([21,20,19,18,17,16,15]):
        ax = plt.subplot(subs[k])
        species = species_in[time,:,:]
        u = u_in[time,:,:]
        v = v_in[time,:,:]
        m.contourf(lon,lat,species,vmin=species_min,vmax=species_max,levels=np.linspace(np.min(species,axis=None),np.max(species,axis=None),301),latlon=True,cmap=matplotlib.cm.jet)
        #m.pcolormesh(lon,lat,species,vmin=species_min,vmax=species_max,latlon=True,cmap=matplotlib.cm.jet,shading='gouraud')
        if k in [3,6]:
            cb = m.colorbar(format='%.3f')
            cb.ax.tick_params(labelsize=8)
            cb.set_clim(species_min,species_max)
        #cb.set_ticklabels(fmt)
        #ax.set_title("{:d}:00 hrs".format(time))
        plt.annotate("{:d}:00".format(time), xy=(0.02, 0.02), xycoords='axes fraction',fontsize=8)
        #    ax.set_title("Time = {:d}".format(-k))
        velquiver = m.quiver(lon[::downsamp,::downsamp],lat[::downsamp,::downsamp],u[::downsamp,::downsamp],v[::downsamp,::downsamp],latlon=True)
        qk = plt.quiverkey(velquiver, 0.9*m.urcrnrx, 1.04*m.urcrnry, 10, '$10$', labelpos='E', coordinates='data', fontproperties={'size': '8'})
        m.drawcoastlines()
        m.drawstates()
        if k in [0,1,2,3]:
            m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=8)
        else:
            m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)
        if k in [0,4]:
            m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
        else:
            m.drawparallels(parallels,labels=[0,0,0,0],fontsize=8)
            
    #plt.tight_layout()
    plt.subplots_adjust(wspace=0.025,hspace=0.15)
    #plt.savefig(s+'quiver_vert.png',dpi=300)
    plt.savefig(s+'quiver_hori.png',dpi=300)
    #plt.savefig(s+'quiver.eps')
    plt.close('all')
























