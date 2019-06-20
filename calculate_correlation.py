from netCDF4 import Dataset
from hdf5storage import loadmat
import numpy as np
list_ = ['H2O','SO2','ANAJ','NO2','O3','ASO4J','ANH4J']
dx = 300
dy = 300

s1_wind = []
for h in range(35):
    print("wind {:02d}".format(h))
    ncfile="ftle/wind_data_{:02d}.nc".format(h)#"ftle_80m.nc"
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    u = vars['eastward_vel'][21,:,:]
    u[u==999999] = np.nan
    v = vars['northward_vel'][21,:,:]
    v[v==999999] = np.nan
    ydim,xdim = u.shape
    s1 = np.empty([ydim,xdim])
    dudy,dudx = np.gradient(u,dy,dx)
    dvdy,dvdx = np.gradient(v,dy,dx)
    for i in range(ydim):
        for j in range(xdim):
            if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
                Utemp = np.array([u[i, j], v[i, j]])
                U = Utemp/np.linalg.norm(Utemp)
                Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
                S = 0.5*(Grad + np.transpose(Grad))
                s1[i,j] = 3600*np.min(np.linalg.eig(S)[0])
            else:
                s1[i,j] = np.nan
    s1_wind.append(s1.ravel())


corr1 = []
for k, species in enumerate(list_):
    print("wind "+species)
    ncfile="ftle/species_data_"+species+".nc"#"ftle_80m.nc"
    root = Dataset(ncfile,'r') #read the data
    vars = root.variables #dictionary, all variables in dataset\
    u = vars['eastward_vel'][21,:,:]
    u[u==999999] = np.nan
    v = vars['northward_vel'][21,:,:]
    v[v==999999] = np.nan
    ydim,xdim = u.shape
    s1 = np.empty([ydim,xdim])
    dudy,dudx = np.gradient(u,dy,dx)
    dvdy,dvdx = np.gradient(v,dy,dx)
    for i in range(ydim):
        for j in range(xdim):
            if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
                Utemp = np.array([u[i, j], v[i, j]])
                U = Utemp/np.linalg.norm(Utemp)
                Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
                S = 0.5*(Grad + np.transpose(Grad))
                s1[i,j] = 3600*np.min(np.linalg.eig(S)[0])
            else:
                s1[i,j] = np.nan
    corr = []
    for h in range(35):
        f_le = s1_wind[h]
        sp_f = s1.ravel()
        
        sp_f = sp_f[~np.isnan(f_le)]
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
    corr1.append(corr)

t=-2
corr2 = []
list_ = ['h2o','so2','anaj','no2','o3','aso4j','anh4j']
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
    corr2.append(corr)
t=0
corr3 = []
list_ = ['h2o','so2','anaj','no2','o3','aso4j','anh4j']
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
    corr3.append(corr)

z=loadmat('flow_map_gridint_3d_q1.mat')['zz'][0,0,:]

np.savez("correlation",c1=corr1,c2=corr2,c3=corr3,z=z)
#"""
