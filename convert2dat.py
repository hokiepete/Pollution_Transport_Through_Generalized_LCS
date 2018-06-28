from netCDF4 import Dataset
import wrf
import numpy as np
height_level = 17 #0.81 eta level
grid_spacing = 12000 #12km

root = Dataset('subset_wrfout_d01_2011-07-01_00_00_00','r')
vars = root.variables
u = vars['U'][:,height_level,:,:]
v = vars['V'][:,height_level,:,:]
lat = vars['XLAT'][:]
lon = vars['XLAT'][:]
root.close()

root = Dataset('subset_wrfout_d01_2011-07-02_00_00_00','r')
vars = root.variables
u = np.concatenate((u,vars['U'][:,height_level,:,:]))
v = np.concatenate((v,vars['V'][:,height_level,:,:]))
lat = np.concatenate((lat,vars['XLAT'][:]))
lon = np.concatenate((lon,vars['XLAT'][:]))
root.close()
u = wrf.destagger(u[:25,:,:],2)
v = wrf.destagger(v[:25,:,:],1)
latitude = lat[:25,:,:]
longitude = lon[:25,:,:]
dim = u.shape
timeout = np.arange(25)/24.0
land = np.zeros(dim,dtype=int)
del lat, lon
print dim
dataset = Dataset('hosiendata.nc', mode='w', format='NETCDF4_CLASSIC') 
lat = dataset.createDimension('lat',dim[1])
lon = dataset.createDimension('lon',dim[2])
time = dataset.createDimension('time', None)

times = dataset.createVariable('time',np.float64,('time',),fill_value=999)
lats = dataset.createVariable('lat',np.float64,('time','lat','lon',),fill_value=999)
lons = dataset.createVariable('lon',np.float64,('time','lat','lon',),fill_value=999)
uo = dataset.createVariable('eastward_vel',np.float64,('time','lat','lon',),fill_value=999)
vo = dataset.createVariable('northward_vel',np.float64,('time','lat','lon',),fill_value=999)
lando = dataset.createVariable('land',np.int,('time','lat','lon',))


lons.standard_name = 'longitude'
lons.units = 'degree_east'
lons.positive = 'east'
lons._CoordinateAxisType = 'Lon'
lons.axis = 'X'
lons.coordsys = 'geographic'
lons[:] = longitude

lats.standard_name = 'latitude'
lats.units = 'degree_north'
lats.positive = 'up'
lats._CoordinateAxisType = 'Lat'
lats.axis = 'Y'
lats.coordsys = 'geographic'
lats[:] = latitude

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
uo[:] = u

vo.standard_name = 'surface_northward_sea_water_velocity'
vo.long_name = 'surface_northward_sea_water_velocity'
vo.units = 'meter second-1'
vo.coordsys = 'geographic'
vo.positive = 'toward north'
vo.coordinates = 'Longitude Latitude datetime'
vo[:] = v

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