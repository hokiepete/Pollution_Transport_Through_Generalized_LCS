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
#label = ['H2O', 'SO2', 'O3', 'SO4', 'NO2', 'Na']
label = ['1500','1800','2100']
subs = [131,132,133]
f = np.load('ridges_flow_Q.npz')
a_r = f['attracting_ridges']
r_r = f['repelling_ridges']
f.close()
for a, s in enumerate(['Q']):#,'SO2','O3','ASO4J','NO2','ANAJ']):
    min_l = []
    max_l = []
    for time in enumerate([15]):#,18,21]):
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
    fig = plt.figure()#1,figsize=FigSize)
    for k, time in enumerate([15,18,21]):
        #print(k,s,subs[k])
        plt.subplot(subs[k])
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
        x,y = m(lon,lat)#lonlat2m(lon=lon,lat=lat,ref_lon=ref_lon,ref_lat=ref_lat)
        
        test = os.listdir()

        for item in test:
            if item.endswith(".png"):
                os.remove(os.path.join(item))
        
        m.pcolormesh(x,y,species,vmin=c_min,vmax=c_max,cmap='Greys',shading='gouraud')
        m.drawcoastlines()
        m.drawstates()
        m.colorbar()
        flow_y = r_r[k]
        flow_x = r_r[k+3]
        for p in range(len(flow_x)):
            x = flow_x[p]
            y = flow_y[p]
            if len(x) > 6:
                m.plot(x,y,'r-')#, latlon=True)
                
        flow_y = a_r[k]
        flow_x = a_r[k+3]
        for p in range(len(flow_x)):
            x = flow_x[p]
            y = flow_y[p]
            if len(x) > 6:
                m.plot(x,y,'b-')#, latlon=True)
        
        
        























