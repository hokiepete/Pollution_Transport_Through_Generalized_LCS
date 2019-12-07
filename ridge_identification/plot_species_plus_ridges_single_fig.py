import netCDF4 as nc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
from mapping_functions import lonlat2m
os.environ["PROJ_LIB"] = "C:\\ProgramData\\Anaconda3\\Library\\share"
from mpl_toolkits.basemap import Basemap

plt.close('all')
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['lines.linewidth'] = 1
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
FigSize=(figwidth, 0.35*ydim/xdim*figwidth)
#.525*
subs = [131,132,133]# ,234,235,236]
#subs = [141,142,143]# ,234,235,236]
#label = ['H2O', 'SO2', 'O3', 'SO4', 'NO2', 'Na']
label = ['1500','1800','2100']
ttimes = [15,18,21]
percen = 85
for a, s in enumerate(['O3','Q','ASO4J','ANAJ']):
    f = np.load(f'ridges_flow_{s}.npz')
    a_r = f['attracting_ridges']
    r_r = f['repelling_ridges']
    f.close()
    min_l = []
    max_l = []
    for time in enumerate([15,18,21]):
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
        max_l.append(species.max())
        min_l.append(species.min())
    
    c_min = min(min_l)
    c_max = max(max_l)
    fig = plt.figure(1,figsize=FigSize)
    for k, time in enumerate(ttimes):
        #print(k,s,subs[k])
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
        
        downsamp = 25
        ref_lat = lat.mean()
        ref_lon = lon.mean()
        lon,lat = np.meshgrid(lon,lat)
        x,y = lonlat2m(lon=lon,lat=lat,ref_lon=ref_lon,ref_lat=ref_lat)
        
        im = m.pcolormesh(lon,lat,species,vmin=c_min,vmax=c_max,latlon=True,cmap='Greys',shading='gouraud')

        flow_y = r_r[k]
        flow_x = r_r[k+3]
        for p in range(len(flow_x)):
            x = flow_x[p]
            y = flow_y[p]
            if len(x) > 7:
                m.plot(x,y,'r-')#, latlon=True)
                
        flow_y = a_r[k]
        flow_x = a_r[k+3]
        for p in range(len(flow_x)):
            x = flow_x[p]
            y = flow_y[p]
            if len(x) > 7:
                m.plot(x,y,'b-')#, latlon=True)

        plt.annotate('{0}'.format(label[k]),xy=(0.03*m.urcrnrx,0.03*m.urcrnry),fontsize=7)
        m.drawcoastlines()
        m.drawstates()
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=7)
        if k ==0:
            m.drawparallels(parallels,labels=[1,0,0,0],fontsize=7)
        else:
            m.drawparallels(parallels,labels=[0,0,0,0],fontsize=7)
    
    cbar_ax = fig.add_axes([0.908, 0.117, 0.01, 0.756])
    cbar = fig.colorbar(im, cax=cbar_ax)
            
    cbar.ax.tick_params(labelsize=6) 
    tick_locator = ticker.MaxNLocator(nbins=9)
    cbar.locator = tick_locator
    cbar.update_ticks()
    plt.subplots_adjust(wspace=0.05, hspace=0)
    plt.savefig(f'{s}_ridges.png',dpi=300)
    plt.close('all')








"""

with np.load(f'{s}_{time}_att_rep.npz') as file:
            s1_directional_derivative=file['arr_0']
            s2_directional_derivative=file['arr_1']
            s1_concavity=file['arr_2']
            s2_concavity=file['arr_3']
            s1=file['arr_4']
            s2=file['arr_5']
        
        s1_directional_derivative = np.ma.masked_where(s1>=0,s1_directional_derivative)
        s2_directional_derivative = np.ma.masked_where(s2<=0,s2_directional_derivative)
        s1_directional_derivative = np.ma.masked_where(s1>=-np.percentile(-s1,percen),s1_directional_derivative)
        s2_directional_derivative = np.ma.masked_where(s2<=np.percentile(s2,percen),s2_directional_derivative)
        
        #Now we apply the third condition, i.e. concavity 
        s1_directional_derivative = np.ma.masked_where(s1_concavity<=0,s1_directional_derivative)
        s2_directional_derivative = np.ma.masked_where(s2_concavity>=0,s2_directional_derivative)
        print(lon.shape)
        im = m.pcolormesh(lon,lat,species,vmin=c_min,vmax=c_max,latlon=True,cmap=matplotlib.cm.jet,shading='gouraud')
        plt.annotate('{0}'.format(label[k]),xy=(0.03*m.urcrnrx,0.03*m.urcrnry),fontsize=7)

        m.contour(lon,lat,s2_directional_derivative,latlon=True,levels=[0],colors='red',linewidths=1.0)
        m.contour(lon,lat,s1_directional_derivative,latlon=True,levels=[0],colors='blue',linewidths=1.0)
        

"""











