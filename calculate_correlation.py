from netCDF4 import Dataset
from hdf5storage import loadmat
import numpy as np
import matplotlib.pyplot as plt
t=-2
list_ = ['03','h2o','so2','anaj','no2','o3','aso4j','anh4j']
z=loadmat('flow_map_gridint_3d_q1.mat')['zz'][0,0,:]




for k, species in enumerate(list_):
    ncfile=species+"_FTLE.nc"#"ftle_80m.nc"
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    time = vars['time'][:]
    sp_ftle = vars['FTLE'][t,:,:].data.ravel()
    sp_ftle[sp_ftle==999999] = np.nan

    corr = []
    for h in range(35):
        ncfile="ftle/{:02d}_ftle.nc".format(h)#"ftle_80m.nc"
        root = Dataset(ncfile,'r') #read the data
        vars = root.variables #dictionary, all variables in dataset\
        f_le = vars['FTLE'][t,:,:].data.ravel()
        f_le[f_le==999999] = np.nan
        
        sp_f = sp_ftle[~np.isnan(f_le)]
        f_le = f_le[~np.isnan(f_le)]
        
        f_le = f_le[~np.isnan(sp_f)]
        sp_f = sp_f[~np.isnan(sp_f)]
        
        bar_le = np.mean(f_le)
        bar_sp = np.mean(sp_f)
        n=len(f_le)
        
        numerator = sum(f_le*sp_f)-(n*bar_le*bar_sp)
        den1 = np.sqrt(sum(f_le**2)-n*bar_le**2)
        den2 = np.sqrt(sum(sp_f**2)-n*bar_sp**2)
        denominator = den1*den2;
        corr.append(numerator/denominator)










for k, species in enumerate(list_):
    ncfile='ftle/'+species+"_FTLE.nc"#"ftle_80m.nc"
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    time = vars['time'][:]
    sp_ftle = vars['FTLE'][t,:,:].data.ravel()
    sp_ftle[sp_ftle==999999] = np.nan

    corr = []
    for h in range(35):
        ncfile="ftle/{:02d}_ftle.nc".format(h)#"ftle_80m.nc"
        root = Dataset(ncfile,'r') #read the data
        vars = root.variables #dictionary, all variables in dataset\
        f_le = vars['FTLE'][t,:,:].data.ravel()
        f_le[f_le==999999] = np.nan
        
        sp_f = sp_ftle[~np.isnan(f_le)]
        f_le = f_le[~np.isnan(f_le)]
        
        f_le = f_le[~np.isnan(sp_f)]
        sp_f = sp_f[~np.isnan(sp_f)]
        
        bar_le = np.mean(f_le)
        bar_sp = np.mean(sp_f)
        n=len(f_le)
        
        numerator = sum(f_le*sp_f)-(n*bar_le*bar_sp)
        den1 = np.sqrt(sum(f_le**2)-n*bar_le**2)
        den2 = np.sqrt(sum(sp_f**2)-n*bar_sp**2)
        denominator = den1*den2;
        corr.append(numerator/denominator)
               
    plt.close('all')
    plt.figure()
    x=np.linspace(0,34,35)
    plt.plot(corr,z)
    plt.savefig('{}.png'.format(species))



"""
from netCDF4 import Dataset
from hdf5storage import loadmat
import numpy as np
f = loadmat('FTLE_wrf_3d.mat')# as f:
FTLE_3d = f['sigma'][:,:,:,6]
del f

    
f= loadmat('FTLE_wrf_3d_levels.mat')# as f:
FTLE = f['sigma'][:,:,:,6]
del f

ftle_species = ['wind_speed','water_vapor','NO2','O3']

ncfile=ftle_species[1]+"_FTLE.nc"#"ftle_80m.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
sp_ftle = vars['FTLE'][-7,:,:,0].data.ravel()
sp_ftle[sp_ftle==999999] = np.nan

corr0 = []
for h in range(35):
    f_3d = FTLE_3d[:,:,h].ravel()
    f_le = FTLE[:,:,h].ravel()
    
    sp_f = sp_ftle[~np.isnan(f_3d)]
    f_3d = f_3d[~np.isnan(f_3d)]
    
    f_3d = f_3d[~np.isnan(sp_f)]
    sp_f = sp_f[~np.isnan(sp_f)]

    bar_3d = np.mean(f_3d)
    bar_sp = np.mean(sp_f)
    n=len(f_3d)
    
    numerator = sum(sp_f*f_3d)-(n*bar_sp*bar_3d)
    den1 = np.sqrt(sum(sp_f**2)-n*bar_sp**2)
    den2 = np.sqrt(sum(f_3d**2)-n*bar_3d**2)
    denominator = den1*den2;
    corr0.append(numerator/denominator)
    
    
    sp_f = sp_ftle[~np.isnan(f_le)]
    f_le = f_le[~np.isnan(f_le)]
    
    f_le = f_le[~np.isnan(sp_f)]
    sp_f = sp_f[~np.isnan(sp_f)]
    
    bar_le = np.mean(f_le)
    bar_sp = np.mean(sp_f)
    n=len(f_le)
    
    numerator = sum(f_le*sp_f)-(n*bar_le*bar_sp)
    den1 = np.sqrt(sum(f_le**2)-n*bar_le**2)
    den2 = np.sqrt(sum(sp_f**2)-n*bar_sp**2)
    denominator = den1*den2;
    corr1.append(numerator/denominator)
    
del FTLE, FTLE_3d

import matplotlib.pyplot as plt

plt.plot(corr)


"""