import mapping_functions as mf
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import calendar
import time
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interp1d
plt.close('all')


height_level = 3 #roughly 80 m above ground level
grid_spacing = 12*1000 #km 2 m

tdim = 25
xdim = 102
ydim = 82

root = Dataset('wrf_2011_07_01','r')
cen_lat = getattr(root,'CEN_LAT')
cen_lon = getattr(root,'CEN_LON')
true_lat1 = getattr(root,'TRUELAT1')
true_lat2 = getattr(root,'TRUELAT2')
ref_lat = getattr(root,'MOAD_CEN_LAT')
ref_lon = getattr(root,'STAND_LON')
vars = root.variables
u = vars['U'][:,height_level,:,:]
v = vars['V'][:,height_level,:,:]
lat_vel = vars['XLAT'][0,:,:]
lon_vel = vars['XLONG'][0,:,:]
root.close()

checklon, checklat = mf.lonlat2km(ref_lon,ref_lat,lon_vel,lat_vel,true_lat1,true_lat2) 

root = Dataset('wrf_2011_07_02','r')
vars = root.variables
u = np.concatenate((u,vars['U'][:,height_level,:,:]))
v = np.concatenate((v,vars['V'][:,height_level,:,:]))
root.close()
u = mf.unstagger(u[:tdim,:,:],2)
v = mf.unstagger(v[:tdim,:,:],1)
dx = grid_spacing
dy = grid_spacing
time_in = np.linspace(0,24,25)
time_want = np.linspace(0,24,145)

u_out = np.empty([145,82,102])
v_out = np.empty([145,82,102])
for i in range(ydim):
    for j in range(xdim):
        fu = interp1d(time_in,u[:,i,j],kind='cubic')
        u_out[:,i,j] = fu(time_want)
        fv = interp1d(time_in,v[:,i,j],kind='cubic')
        v_out[:,i,j] = fv(time_want)

u = u_out
v = v_out
del u_out, v_out

s1 = np.ma.empty(u.shape)

J = np.array([[0, 1], [-1, 0]])
for t in range(u.shape[0]):
    dudy,dudx = np.gradient(u[t,:,:],dy,dx,edge_order=2)
    dvdy,dvdx = np.gradient(v[t,:,:],dy,dx,edge_order=2)
    for i in range(ydim):
        for j in range(xdim):
            if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[t,i,j] and v[t,i,j]) is not np.ma.masked:    
                Utemp = np.array([u[t,i, j], v[t,i, j]])
                Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
                S = 0.5*(Grad + np.transpose(Grad))
                eigenValues, eigenVectors_temp = np.linalg.eig(S)
                idx = eigenValues.argsort()
                s1[t,i,j] = eigenValues[idx[1]]
    
            else:
                s1[t,i,j] = np.ma.masked



xdim = 405
ydim = 325
tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
ncfile="ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
lat = vars["lat"][:]#.reshape([ydim,xdim])
lon = vars["lon"][:]#.reshape([ydim,xdim])
print(vars['FTLE'][:].shape)
time = 86400*vars["time"][:]+tstart
tdim = np.shape(time)[0]
#ftle = vars['FTLE'][:,:,:,0]#.reshape([ydim,xdim,tdim],order='C')
root.close()
lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)

m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'c',
            area_thresh=1000.,
            )
'''
m = Basemap(llcrnrlon=-179,
            llcrnrlat=-85,
            urcrnrlon=180,
            urcrnrlat=85,
            projection='merc',
            resolution = 'c',
            area_thresh=1000.,
            )
#'''


#lon,lat = np.meshgrid(lon,lat)
#cs=m.contourf(lon,lat,ftle[-1,:,:],levels=np.linspace(ftle.min(axis=None),ftle.max(axis=None),301),latlon=True)
#ncfile="SE_tracers.nc"
ncfile="tracers.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
t_lat = vars["lat"][:]#.reshape([ydim,xdim])
t_lon = vars["lon"][:]#.reshape([ydim,xdim])
time = 86400*vars["time"][:]+tstart
tdim = np.shape(time)[0]
root.close()
attracting_ridges_lat=[]
attracting_ridges_lon=[]
for a in [13,17,40,46,47,55]:
    ncfile="attracting_ridge_{:02d}.nc".format(a)
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    attracting_ridges_lat.append(vars["lat"][:])
    attracting_ridges_lon.append(vars["lon"][:])
    root.close()

repelling_ridges_lat=[]
repelling_ridges_lon=[]
for r in [8,18,19,53,56,60]:
    ncfile="repelling_ridge_{:02d}.nc".format(r)
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    repelling_ridges_lat.append(vars["lat"][:])
    repelling_ridges_lon.append(vars["lon"][:])
    root.close()
    
for t in range(time.shape[0]):
    #for c in cs.collections:
    #    c.remove()
    #cs.set_array(np.ravel(ftle[:,:,t]))
    #cs=m.contourf(lon,lat,ftle[:,:,-t],levels=np.linspace(ftle.min(axis=None),ftle.max(axis=None),301),latlon=True)

    #cs=m.contourf(lon,lat,ftle[-t,:,:],levels=np.linspace(ftle[-t,:,:].min(axis=None),ftle[-t,:,:].max(axis=None),301),latlon=True)
    #plt.title("{0}".format(time[-t]),fontsize=18)
    hrs, mins = np.divmod(t*10,60)
    plt.figure(figsize=[16,12])
    m.contourf(lon_vel,lat_vel,s1[t,:,:],levels=np.linspace(s1[t,:,:].min(axis=None),s1[t,:,:].max(axis=None),301),latlon=True)
    #m.contourf(lon_vel,lat_vel,-s1[t,:,:],levels=np.linspace(-s1[t,:,:].max(axis=None),-s1[t,:,:].min(axis=None),301),latlon=True)
    for i in range(len(attracting_ridges_lat)):
        x = [elem for elem in attracting_ridges_lon[i][:,t].data if elem < 99999]
        y = [elem for elem in attracting_ridges_lat[i][:,t].data if elem < 99999]
        #m.plot(attracting_ridges_lon[i][:,t],attracting_ridges_lat[i][:,t],c='b',latlon=True)
        m.plot(x,y,c='b',latlon=True)
    for i in range(len(repelling_ridges_lat)):
        x = [elem for elem in repelling_ridges_lon[i][:,t].data if elem < 99999]
        y = [elem for elem in repelling_ridges_lat[i][:,t].data if elem < 99999]
        #m.plot(repelling_ridges_lon[i][:,t],repelling_ridges_lat[i][:,t],c='r',latlon=True)
        m.plot(x,y,c='r',latlon=True)
    m.scatter(t_lon[:,t],t_lat[:,t],latlon=True)
    plt.title('{0:02d} hrs, {1:02d} mins'.format(int(hrs),int(mins)))
    
    m.drawcoastlines()
    m.drawstates()
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-2)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)    
    plt.savefig('OECS+s2_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
    
    plt.close('all')
#'''
    
import winsound
frequency = 800  # Set Frequency To 2500 Hertz
duration = 500  # Set Duration To 1000 ms == 1 second
winsound.Beep(frequency, duration)
    
"""
ncfile="SE_ridge.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
r_lat = vars["lat"][:]#.reshape([ydim,xdim])
r_lon = vars["lon"][:]#.reshape([ydim,xdim])
root.close()

"""