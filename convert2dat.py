from netCDF4 import Dataset
import wrf
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
    #stdlat1  =  deg2rad(30)
    #stdlat2  =  deg2rad(60)
    stdlat1  =  deg2rad(40)
    stdlat2  =  deg2rad(42)
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
u = vars['U'][:,height_level,:,:]
v = vars['V'][:,height_level,:,:]
lat = vars['XLAT'][0,:,:]
lon = vars['XLONG'][0,:,:]
root.close()

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

xout,yout = lonlat2km(0.5*(np.min(lon)+np.max(lon)),0.5*(np.min(lat)+np.max(lat)),lonout,latout)
uout = np.empty([tdim,yodim*xodim])
vout = np.empty([tdim,yodim*xodim])

for t in range(tdim):
    uout[t,:] = interpolate.griddata((yin.ravel(), xin.ravel()), u[t,:,:].ravel(),(yout.ravel(), xout.ravel()),method='cubic',fill_value=999)
    vout[t,:] = interpolate.griddata((yin.ravel(), xin.ravel()), v[t,:,:].ravel(),(yout.ravel(), xout.ravel()),method='cubic',fill_value=999)

uout = np.reshape(uout, [tdim,yodim,xodim])
vout = np.reshape(vout, [tdim,yodim,xodim])


dim = uout.shape
timeout = np.arange(tdim)/24.0
land = np.zeros(dim,dtype=int)

print dim
dataset = Dataset('hosiendata.nc', mode='w', format='NETCDF4_CLASSIC') 
lat = dataset.createDimension('lat',dim[1])
lon = dataset.createDimension('lon',dim[2])
time = dataset.createDimension('time', None)

times = dataset.createVariable('time',np.float64,('time',),fill_value=999)
lats = dataset.createVariable('lat',np.float64,('lat'),fill_value=999)
lons = dataset.createVariable('lon',np.float64,('lon',),fill_value=999)
uo = dataset.createVariable('eastward_vel',np.float64,('time','lat','lon',),fill_value=999)
vo = dataset.createVariable('northward_vel',np.float64,('time','lat','lon',),fill_value=999)
lando = dataset.createVariable('land',np.int,('time','lat','lon',))


lons.standard_name = 'longitude'
lons.units = 'degree_east'
lons.positive = 'east'
lons._CoordinateAxisType = 'Lon'
lons.axis = 'X'
lons.coordsys = 'geographic'
lons[:] = lon2file

lats.standard_name = 'latitude'
lats.units = 'degree_north'
lats.positive = 'up'
lats._CoordinateAxisType = 'Lat'
lats.axis = 'Y'
lats.coordsys = 'geographic'
lats[:] = lat2file

times.standard_name = 'time'
times.long_name = 'time'
times.units = 'days since 2011-07-01 00:00:00 UTC'
times.calendar = 'gregorian'
times._CoordinateAxisType = 'Time'
times[:] = timeout[:]

uo.standard_name = 'surface_eastward_sea_water_velocity'
uo.long_name = 'surface_eastward_sea_water_velocity'
uo.units = 'meter second-1'
uo.coordsys = 'geographic'
uo.positive = 'toward east'
uo.coordinates = 'Longitude Latitude datetime'
uo[:] = uout

vo.standard_name = 'surface_northward_sea_water_velocity'
vo.long_name = 'surface_northward_sea_water_velocity'
vo.units = 'meter second-1'
vo.coordsys = 'geographic'
vo.positive = 'toward north'
vo.coordinates = 'Longitude Latitude datetime'
vo[:] = vout

lando.standard_name = 'land'
lando.units = 'Boolean'
lando.coordinates = 'Longitude Latitude datetime'
lando[:] = land

dataset.close()


'''
dim = u.shape
xmin = -0.5*(dim[2]-1)*grid_spacing
xmax = 0.5*(dim[2]-1)*grid_spacing
ymin = -0.5*(dim[1]-1)*grid_spacing
ymax = 0.5*(dim[1]-1)*grid_spacing
umax = n.max(n.abs(u))*3600*24 - (dim[2]-1)*grid_spacing
vmax = n.max(n.abs(v))*3600*24 - (dim[1]-1)*grid_spacing
'''