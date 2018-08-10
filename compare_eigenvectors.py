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
            s1[i,j] = eigenValues[idx[0]]
            eigenMax[i,j,:] = eigenVectors_temp[:,idx[-1]]
            s2[i,j] = eigenValues[idx[-1]]

        else:
            s1[i,j] = np.ma.masked
            eigenMin[i,j,0] = np.ma.masked
            eigenMin[i,j,1] = np.ma.masked
            s2[i,j] = np.ma.masked
            eigenMax[i,j,0] = np.ma.masked
            eigenMax[i,j,1] = np.ma.masked
