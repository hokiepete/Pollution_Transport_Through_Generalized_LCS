from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
from scipy import interpolate
def cot(th):
    return 1.0/np.tan(th)

def sec(th):
    return 1.0/np.cos(th)

def deg2rad(deg):
    return deg*np.pi/180.0

def rad2deg(rad):
    return rad*180.0/np.pi

def lonlat2km(reflon,reflat,lon,lat ):
    #LONLAT2KM Summary of this function goes here
    #   Uses Lambert Conformal Projection
    stdlat1  =  deg2rad(30)
    stdlat2  =  deg2rad(60)
    R=6371
    reflon = deg2rad(reflon)
    reflat = deg2rad(reflat)
    lon = deg2rad(lon)
    lat = deg2rad(lat)
    n = np.log(np.cos(stdlat1)*sec(stdlat2)) / np.log(np.tan(0.25*np.pi+0.5*stdlat2)*cot(0.25*np.pi+0.5*stdlat1))
    F=(np.cos(stdlat1)*(np.tan(0.25*np.pi+0.5*stdlat1)**n))/n
    p0 = R*F*(cot(0.25*np.pi+0.5*reflat)**n)
    p = R*F*(cot(0.25*np.pi+0.5*lat)**n)
    th = n*(lon-reflon)
    x=p*np.sin(th)
    y=p0-p*np.cos(th)
    return x,y


height_level = 17 #0.81 eta level
grid_spacing = 12 #km

tdim = 25
xdim = 102
ydim = 82


root = Dataset('subset_wrfout_d01_2011-07-01_00_00_00','r')
vars = root.variables
grd = vars['HGT'][0,:,:]
hgt = (vars['PH'][0,:,:,:] + vars['PHB'][0,:,:,:])/9.81
lat = vars['XLAT'][0,:,:]
lon = vars['XLONG'][0,:,:]
root.close()

import matplotlib.pyplot as plt
data = hgt[3,:,:]-grd
data = np.ma.masked_where(data<0,data)
plt.close('all')
plt.contourf(lon,lat,data,levels=np.linspace(data.min(),data.max(),301))
plt.colorbar()


'''

root = Dataset('subset_wrfout_d01_2011-07-02_00_00_00','r')
vars = root.variables
u = np.concatenate((u,vars['U'][:,height_level,:,:]))
v = np.concatenate((v,vars['V'][:,height_level,:,:]))
latin = np.concatenate((lat,vars['XLAT'][0,:,:]))
lonin = np.concatenate((lon,vars['XLONG'][0,:,:]))
root.close()
u = wrf.destagger(u[:25,:,:],2)
v = wrf.destagger(v[:25,:,:],1)
#latin = lat[:25,:,:]
#longin = lon[:25,:,:]

xin = np.linspace(-grid_spacing*(xdim-1)/2,grid_spacing*(xdim-1)/2,xdim)
yin = np.linspace(-grid_spacing*(ydim-1)/2,grid_spacing*(ydim-1)/2,ydim)
xin, yin = np.meshgrid(xin,yin)


lat2file = np.linspace(np.min(lat),np.max(lat),(ydim-1)*4+1)
lon2file = np.linspace(np.min(lon),np.max(lon),(xdim-1)*4+1)
xodim = lon2file.shape[0]
yodim = lat2file.shape[0]
lonout, latout = np.meshgrid(lon2file,lat2file)

'''