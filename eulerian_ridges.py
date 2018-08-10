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
'''
u=u[-1,:,:]
v=v[-1,:,:]
'''
u=u[0,:,:]
v=v[0,:,:]
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
eigenMin = np.ma.empty([ydim,xdim,2])
s2 = np.ma.empty([ydim,xdim])
eigenMax = np.ma.empty([ydim,xdim,2])
#det = np.ma.empty([ydim,xdim])

J = np.array([[0, 1], [-1, 0]])
for i in range(ydim):
    for j in range(xdim):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Utemp = np.array([u[i, j], v[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            eigenValues, eigenVectors_temp = np.linalg.eig(S)
            idx = eigenValues.argsort()
            eigenMin[i,j,:] = eigenVectors_temp[:,idx[0]]
            #eigenMin[i,j,0] = eigenVectors_temp[1,idx[0]]
            s1[i,j] = eigenValues[idx[0]]
            eigenMax[i,j,:] = eigenVectors_temp[:,idx[-1]]
            #eigenMax[i,j,0] = eigenVectors_temp[1,idx[-1]]
            s2[i,j] = eigenValues[idx[-1]]

        else:
            s1[i,j] = np.ma.masked
            eigenMin[i,j,0] = np.ma.masked
            eigenMin[i,j,1] = np.ma.masked
            s2[i,j] = np.ma.masked
            eigenMax[i,j,0] = np.ma.masked
            eigenMax[i,j,1] = np.ma.masked
            

adfdy,adfdx = np.gradient(s1,dy,dx,edge_order=2)
adfdydy,adfdydx = np.gradient(adfdy,dy,dx,edge_order=2)
adfdxdy,adfdxdx = np.gradient(adfdx,dy,dx,edge_order=2)
adirdiv = np.ma.empty([ydim,xdim])
aconcav = np.ma.empty([ydim,xdim])


rdfdy,rdfdx = np.gradient(s2,dy,dx,edge_order=2)
rdfdydy,rdfdydx = np.gradient(rdfdy,dy,dx,edge_order=2)
rdfdxdy,rdfdxdx = np.gradient(rdfdx,dy,dx,edge_order=2)
rdirdiv = np.ma.empty([ydim,xdim])
rconcav = np.ma.empty([ydim,xdim])

for i in range(ydim):
    for j in range(xdim):
        if (adfdx[i,j] and adfdy[i,j] and adfdxdy[i,j] and adfdydy[i,j] and adfdxdx[i,j] and adfdydx[i,j]) is not np.ma.masked:    
            adirdiv[i,j] = np.dot([adfdx[i,j],adfdy[i,j]],eigenMin[i,j,:])
            aconcav[i,j] = np.dot(np.dot([[adfdxdx[i,j],adfdxdy[i,j]],[adfdydx[i,j],adfdydy[i,j]]],eigenMin[i,j,:]),eigenMin[i,j,:])
        #print aconcav[i,j]
        else:
            adirdiv[i,j] = np.ma.masked
            aconcav[i,j] = np.ma.masked

        if (rdfdx[i,j] and rdfdy[i,j] and rdfdxdy[i,j] and rdfdydy[i,j] and rdfdxdx[i,j] and rdfdydx[i,j]) is not np.ma.masked:    
            rdirdiv[i,j] = np.dot([rdfdx[i,j],rdfdy[i,j]],eigenMax[i,j,:])
            rconcav[i,j] = np.dot(np.dot([[rdfdxdx[i,j],rdfdxdy[i,j]],[rdfdydx[i,j],rdfdydy[i,j]]],eigenMax[i,j,:]),eigenMax[i,j,:])
        #print aconcav[i,j]
        else:
            rdirdiv[i,j] = np.ma.masked
            rconcav[i,j] = np.ma.masked


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

m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'i',
            area_thresh=1000.,
            )
#'''
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)

fig = plt.figure(1)
plt.subplot(121)
#lon,lat = np.meshgrid(lon,lat)
cs = m.contourf(lon,lat,s1,levels=np.linspace(np.min(s1,axis=None),np.max(s1,axis=None),301),latlon=True)
plt.colorbar()
adirdiv = np.ma.masked_where(aconcav<0,adirdiv)
adirdiv = np.ma.masked_where(-s1<=np.percentile(-s1,85),adirdiv)
ridge = m.contour(lon,lat,adirdiv,levels=0,latlon=True,colors='cyan')
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

plt.subplot(122)
#lon,lat = np.meshgrid(lon,lat)
cs = m.contourf(lon,lat,s2,levels=np.linspace(np.min(s2,axis=None),np.max(s2,axis=None),301),latlon=True)
plt.colorbar()
rdirdiv = np.ma.masked_where(rconcav>0,rdirdiv)
rdirdiv = np.ma.masked_where(s2<=np.percentile(s2,85),rdirdiv)
ridge = m.contour(lon,lat,rdirdiv,levels=0,latlon=True,colors='magenta')
m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

fig = plt.figure(2)
m.quiver(lon[::5,::5],lat[::5,::5],u[::5,::5],v[::5,::5],latlon=True)
ridge = m.contour(lon,lat,adirdiv,levels=0,latlon=True,colors='blue')
ridge = m.contour(lon,lat,rdirdiv,levels=0,latlon=True,colors='red')
#ridge = m.contour(lon,lat,det,levels=0,latlon=True)

m.drawcoastlines()
m.drawstates()
parallels = np.arange(round(lat_min,0),lat_max+2,2)
meridians = np.arange(round(lon_max,0),lon_min-2,-2)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

'''
t=1
hrs, mins = np.divmod((t-1)*10,60)
plt.title("Integration time = {0:02d} hrs, {1:02d} min".format(hrs,mins),fontsize=18)
plt.savefig('SE_lcs_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
'''
