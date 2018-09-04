from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
from scipy import interpolate
import mapping_functions as mf
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker as ticker
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}


height_level = 17 #0.81 eta level
height_level = 3 #roughly 80 m above ground level
grid_spacing = 12*1000 #km 2 m

tdim = 25
xdim = 102
ydim = 82

figwidth = 6
FigSize=(figwidth, ydim/xdim*figwidth)


root = Dataset('subset_wrfout_d01_2011-07-01_00_00_00','r')
#root = Dataset('wrf_2011_07_01','r')
#root = Dataset('subset_wrfout_d01_2011-07-02_00_00_00')
cen_lat = getattr(root,'CEN_LAT')
cen_lon = getattr(root,'CEN_LON')
true_lat1 = getattr(root,'TRUELAT1')
true_lat2 = getattr(root,'TRUELAT2')
ref_lat = getattr(root,'MOAD_CEN_LAT')
ref_lon = getattr(root,'STAND_LON')
vars = root.variables
#Wind Velocity
u = vars['U'][:,height_level,:,:]
v = vars['V'][:,height_level,:,:]
#Water Vapor Flux, Vertically Integrated
#u = vars['UQ'][:,:,:]
#v = vars['VQ'][:,:,:]
lat = vars['XLAT'][0,:,:]
lon = vars['XLONG'][0,:,:]
root.close()
checklon, checklat = mf.lonlat2km(ref_lon,ref_lat,lon,lat,true_lat1,true_lat2) 

root = Dataset('subset_wrfout_d01_2011-07-02_00_00_00','r')
vars = root.variables
#Wind Velocity
u = np.concatenate((u,vars['U'][:,height_level,:,:]))
v = np.concatenate((v,vars['V'][:,height_level,:,:]))
#Water Vapor Flux, Vertically Integrated
#u = np.concatenate((u,vars['UQ'][:,:,:]))
#v = np.concatenate((v,vars['VQ'][:,:,:]))
root.close()
u = mf.unstagger(u[:25,:,:],2)
v = mf.unstagger(v[:25,:,:],1)

u=u[-2,:,:]
v=v[-2,:,:]
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
s2 = np.ma.empty([ydim,xdim])
J = np.array([[0, 1], [-1, 0]])
for i in range(ydim):
    for j in range(xdim):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Utemp = np.array([u[i, j], v[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            s1[i,j] = 3600*np.min(np.linalg.eig(S)[0])
            s2[i,j] = 3600*np.max(np.linalg.eig(S)[0])

        else:
            s1[i,j] = np.ma.masked
            s2[i,j] = np.ma.masked

plt.close('all')
lev= 301
parallels_spacing = 3
meridian_spacing = -3
lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)
height=(ydim-1)*dy#+10*1000
width=(xdim-1)*dx#+10*1000

m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'l',
            area_thresh=5000.,
            )
#'''
parallels = np.arange(round(lat_min,0),lat_max+2,parallels_spacing)
meridians = np.arange(round(lon_max,0),lon_min-2,meridian_spacing)

fig = plt.figure(figsize=FigSize)
sub = plt.subplot(221)
cs = m.contourf(lon,lat,-s1,levels=np.linspace(np.min(-s1,axis=None),np.max(-s1,axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
#plt.title("s$_{1}$ field",**titlefont)#fontsize=18)
plt.annotate('A', xy=(0.92, 0.03), xycoords='axes fraction')

ncfile="ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
lat = vars["lat"][:]#.reshape([ydim,xdim])
lon = vars["lon"][:]#.reshape([ydim,xdim])
print(vars['FTLE'][:].shape)
ftle = vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
#lon, lat = np.meshgrid(lon,lat)
root.close()
lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)
lon, lat = np.meshgrid(lon,lat)
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'l',
            area_thresh=5000.,
            )

parallels = np.arange(round(lat_min,0),lat_max+2,parallels_spacing)
meridians = np.arange(round(lon_max,0),lon_min-2,meridian_spacing)
sub = plt.subplot(222)
t=4
cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
#plt.ylabel('hr$^{-1}$',**labelfont)
hrs, mins = np.divmod((t-1)*10,60)
#plt.title("Integration time = -{0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
#plt.title("Integration time = -{0:02d}.{1:02g} hrs".format(hrs,mins/60),**titlefont)
plt.annotate('B', xy=(0.92, 0.03), xycoords='axes fraction')

sub = plt.subplot(223)
t=7
cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
#plt.ylabel('hr$^{-1}$',**labelfont)
hrs, mins = np.divmod((t-1)*10,60)
#plt.title("Integration time = -{0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
#plt.title("Integration time = -{0:02d}.{1:02g} hrs".format(hrs,mins/60),**titlefont)
plt.annotate('C', xy=(0.92, 0.03), xycoords='axes fraction')

sub = plt.subplot(224)
t=13
cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),lev),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],**tickfont)
m.drawmeridians(meridians,labels=[0,0,0,1],**tickfont)
plt.yticks(**tickfont)
plt.xticks(**tickfont)
#plt.ylabel('hr$^{-1}$',**labelfont)
hrs, mins = np.divmod((t-1)*10,60)
#plt.title("Integration time = -{0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
plt.annotate('D', xy=(0.92, 0.03), xycoords='axes fraction')

plt.savefig('s1_Backward-Time_FTLE_Comparison.eps', transparent=False, bbox_inches='tight')
