from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
from scipy import interpolate
import mapping_functions as mf
    
height_level = 17 #0.81 eta level
height_level = 3 #roughly 80 m above ground level
grid_spacing = 12*1000 #km 2 m

tdim = 25
xdim = 102
ydim = 82

root = Dataset('subset_wrfout_d01_2011-07-01_00_00_00','r')
cen_lat = getattr(root,'CEN_LAT')
cen_lon = getattr(root,'CEN_LON')
true_lat1 = getattr(root,'TRUELAT1')
true_lat2 = getattr(root,'TRUELAT2')
ref_lat = getattr(root,'MOAD_CEN_LAT')
ref_lon = getattr(root,'STAND_LON')
vars = root.variables
u = vars['U'][:,height_level,:,:]
v = vars['V'][:,height_level,:,:]
lat = vars['XLAT'][0,:,:]
lon = vars['XLONG'][0,:,:]
root.close()

checklon, checklat = mf.lonlat2km(ref_lon,ref_lat,lon,lat,true_lat1,true_lat2) 

root = Dataset('subset_wrfout_d01_2011-07-02_00_00_00','r')
vars = root.variables
u = np.concatenate((u,vars['U'][:,height_level,:,:]))
v = np.concatenate((v,vars['V'][:,height_level,:,:]))
root.close()
u = mf.unstagger(u[:25,:,:],2)
v = mf.unstagger(v[:25,:,:],1)
u=u[-1,:,:]
v=v[-1,:,:]
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
rhodot = np.ma.empty([ydim,xdim])
nudot = np.ma.empty([ydim,xdim])
J = np.array([[0, 1], [-1, 0]])
for i in range(ydim):
    for j in range(xdim):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Utemp = np.array([u[i, j], v[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            s1[i,j] = 3600*np.min(np.linalg.eig(S)[0])
            s2[i,j] = 3600*np.max(np.linalg.eig(S)[0])
            rhodot[i, j] = 3600*np.dot(Utemp, np.dot(np.dot(np.transpose(J), np.dot(S, J)), Utemp))/np.dot(Utemp, Utemp)
            nudot[i, j] = 3600*np.dot(Utemp,np.dot(np.trace(S)*np.identity(2) - 2*S, Utemp))/np.dot(Utemp, Utemp)

        else:
            rhodot[i,j] = np.ma.masked
            nudot[i,j] = np.ma.masked
            s2[i,j] = np.ma.masked

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker as ticker
plt.register_cmap(name='co', data=mf.co())
plt.close('all')

lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)
#reflat=(lat_max + lat_min)/2,
#reflon=(lon_max + lon_min)/2,
height=(ydim-1)*dy#+10*1000
width=(xdim-1)*dx#+10*1000
'''
m = Basemap(height=height,
            width=width,
            lat_1 = true_lat1,
            lat_2 = true_lat2,
            lat_0=cen_lat,
            lon_0=cen_lon,
            lon_1=ref_lon,
            projection='lcc',
            resolution = 'h',
            area_thresh=1000.,
            )
'''
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'h',
            area_thresh=1000.,
            )
#'''
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)
rholevel = np.min(np.abs([rhodot.min(axis=None),rhodot.max(axis=None)]))
nulevel = np.min(np.abs([nudot.min(axis=None),nudot.max(axis=None)]))
'''
#plt.rcParams(usetex)
fig = plt.figure(1,figsize=(16,12))
ax = plt.subplot(221)
cs = m.contourf(lon,lat,rhodot,levels=np.linspace(rhodot.min(axis=None),rhodot.max(axis=None),301),\
                vmin=-2/3*rholevel,vmax=2/3*rholevel,latlon=True,cmap='co')
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.colorbar()#format=ticker.FuncFormatter(mf.fmt))
plt.title('$\dot{\\rho}$')

ax = plt.subplot(222)
cs = m.contourf(lon,lat,nudot,levels=np.linspace(nudot.min(axis=None),nudot.max(axis=None),301),\
                vmin=-2/3*nulevel,vmax=2/3*nulevel,latlon=True,cmap='PiYG')
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.colorbar()
plt.title('$\dot{\\nu}$')

rdim = rhodot.shape
for i in range(rdim[0]):
    for j in range(rdim[1]):
        #if (rhodot[i,j]>0 and nudot[i,j]>0) or (rhodot[i,j]<0 and nudot[i,j]<0):
        if (rhodot[i,j]<0 and nudot[i,j]<0):
            'nothing'
            #rhodot = rhodot[(rhodot>0 and nudot>0) or (rhodot<0 and nudot<0)]
        else:
            rhodot[i,j] = np.ma.masked

ax = plt.subplot(223)
#cs = m.contourf(lon,lat,-rhodot,levels=np.linspace(np.min(-rhodot,axis=None),np.max(-rhodot,axis=None),301),\
#                vmin=-2/3*rholevel,vmax=2/3*rholevel,latlon=True)
cs = m.contourf(lon,lat,-rhodot,levels=np.linspace(np.min(-rhodot,axis=None),np.max(-rhodot,axis=None),301),\
                latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.colorbar()
plt.title('$\dot{\\rho}$ filtered by $\dot{\\nu}$')


ax = plt.subplot(224)
cs = m.contourf(lon,lat,-s1,levels=np.linspace(np.min(-s1,axis=None),np.max(-s1,axis=None),301),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.colorbar()
plt.title('$s_{1}$')

plt.show()

fig = plt.figure(2,figsize=(16,12))
ax = plt.subplot(111)
cs = m.contourf(lon,lat,s2,levels=np.linspace(s2.min(axis=None),s2.max(axis=None),301),latlon=True)
m.drawcoastlines()
m.drawstates()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.colorbar()
plt.title('$s_{2}$')

plt.show()

'''
#fig = plt.figure(3)
#lon,lat = np.meshgrid(lon,lat)
cs = m.contourf(lon,lat,-s1,levels=np.linspace(np.min(-s1,axis=None),np.max(-s1,axis=None),301),latlon=True)
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
t=1
hrs, mins = np.divmod((t-1)*10,60)
plt.title("Integration time = {0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
plt.savefig('SE_lcs_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')

