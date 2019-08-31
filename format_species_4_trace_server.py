from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
from scipy import interpolate
import mapping_functions as mf
    
height_level = 17 #0.81 eta level
height_level = 3 #roughly 80 m above ground level
grid_spacing = 12 #km

tdim = 24
xdim = 102
ydim = 82

for species in ['NO2','SO2','O3','ANH4J','ASO4J','ANAJ']:
    root = Dataset('wrf_species_2011_07_01','r')
    cen_lat = getattr(root,'CEN_LAT')
    cen_lon = getattr(root,'CEN_LON')
    true_lat1 = getattr(root,'TRUELAT1')
    true_lat2 = getattr(root,'TRUELAT2')
    ref_lat = getattr(root,'MOAD_CEN_LAT')
    ref_lon = getattr(root,'STAND_LON')
    vars = root.variables
    #Wind Velocity
    #u = vars['U'][:,height_level,:,:]
    #v = vars['V'][:,height_level,:,:]
    #Water Vapor Flux, Vertically Integrated
    u = vars['U'+species][:,1:-1,1:-1]
    v = vars['V'+species][:,1:-1,1:-1]
    sp = np.mean(vars[species][:,:,1:-1,1:-1],axis=1)   
    lat_in = vars['XLAT'][0,1:-1,1:-1]
    lon_in = vars['XLONG'][0,1:-1,1:-1]
    
    root.close()
    
    xin,yin = mf.lonlat2km(ref_lon,ref_lat,lon_in,lat_in,true_lat1,true_lat2)
    
    lat2file = np.linspace(np.min(lat_in),np.max(lat_in),(ydim-1)*4+1)
    lon2file = np.linspace(np.min(lon_in),np.max(lon_in),(xdim-1)*4+1)
    xodim = lon2file.shape[0]
    yodim = lat2file.shape[0]
    lonout, latout = np.meshgrid(lon2file,lat2file)
    
    xout,yout = mf.lonlat2km(ref_lon,ref_lat,lonout,latout,true_lat1,true_lat2)
    
    uout = np.empty([tdim,yodim*xodim])
    vout = np.empty([tdim,yodim*xodim])
    spout = np.empty([tdim,yodim*xodim])
    
    for t in range(tdim):
        uout[t,:] = interpolate.griddata((yin.ravel(), xin.ravel()), u[t,:,:].ravel(),(yout.ravel(), xout.ravel()),method='cubic',fill_value=999)
        vout[t,:] = interpolate.griddata((yin.ravel(), xin.ravel()), v[t,:,:].ravel(),(yout.ravel(), xout.ravel()),method='cubic',fill_value=999)
        spout[t,:] = interpolate.griddata((yin.ravel(), xin.ravel()), sp[t,:,:].ravel(),(yout.ravel(), xout.ravel()),method='cubic',fill_value=999)
    
    uout = np.reshape(uout, [tdim,yodim,xodim])
    vout = np.reshape(vout, [tdim,yodim,xodim])
    spout = np.reshape(spout, [tdim,yodim,xodim])
    
    dim = uout.shape
    timeout = np.arange(tdim)/23.0
    land = np.zeros(dim,dtype=int)
    
    print(dim)
    dataset = Dataset('species_data_'+species+'.nc', mode='w', format='NETCDF4_CLASSIC') 
    lat = dataset.createDimension('lat',dim[1])
    lon = dataset.createDimension('lon',dim[2])
    time = dataset.createDimension('time', None)
    
    times = dataset.createVariable('time',np.float64,('time',),fill_value=999)
    lats = dataset.createVariable('lat',np.float64,('lat'),fill_value=999)
    lons = dataset.createVariable('lon',np.float64,('lon',),fill_value=999)
    uo = dataset.createVariable('eastward_vel',np.float64,('time','lat','lon',),fill_value=999)
    vo = dataset.createVariable('northward_vel',np.float64,('time','lat','lon',),fill_value=999)
    #This line for summer 2019 only.
    sp = dataset.createVariable('species',np.float64,('time','lat','lon',),fill_value=999)
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
    
    sp.standard_name = 'The species'
    sp.long_name = 'The species'
    sp.units = 'concentration'
    sp.coordsys = 'geographic'
    sp.positive = 'better question would be whats negative mean'
    sp.coordinates = 'Longitude Latitude datetime'
    sp[:] = spout
    
    lando.standard_name = 'land'
    lando.units = 'Boolean'
    lando.coordinates = 'Longitude Latitude datetime'
    lando[:] = land
    
    dataset.close()