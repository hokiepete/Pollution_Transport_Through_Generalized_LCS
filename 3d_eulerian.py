# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 12:07:59 2018

@author: pnola
"""

from netCDF4 import Dataset
import numpy as np
import mapping_functions as f
#from mpl_toolkits.basemap import interp
from scipy import interpolate
root = Dataset('wrf_2011_07_01','r')
vars = root.variables
grd = vars['HGT'][0,:,:]
hgt = (vars['PH'][0,:,:,:] + vars['PHB'][0,:,:,:])/9.81
#hgt = f.unstagger(hgt,axis=0)
lat = vars['XLAT'][0,:,:]
lon = vars['XLONG'][0,:,:]
u = f.unstagger(vars['U'][:],axis=3)
v = f.unstagger(vars['V'][:],axis=2)
w = f.unstagger(vars['W'][:],axis=1)
root.close()
hgt=hgt-grd
heightmax=hgt.max(axis=(1,2))
heightmin=hgt.min(axis=(1,2))
heightmean=hgt.mean(axis=(1,2))

dx=12000
dy=12000
dz = np.diff(heightmean)
xx = np.linspace(0,dx*(hgt.shape[2]-1),hgt.shape[2])
yy = np.linspace(0,dy*(hgt.shape[1]-1),hgt.shape[1])
[yy,zz,xx] = np.meshgrid(yy,heightmean,xx)
t=0
[dudz,dudy,dudx] = np.gradient(u[t,:,:,:],dz,dy,dx)
[dvdz,dvdy,dvdx] = np.gradient(v[t,:,:,:],dz,dy,dx)
[dwdz,dwdy,dwdx] = np.gradient(w[t,:,:,:],dz,dy,dx)
s1 = np.empty(dudz.shape)
for k in range(dudz.shape[0]):
    for i in range(dudz.shape[1]):
        for j in range(dudz.shape[2]):
            dF=np.array([[dudx[k,i,j],dudy[k,i,j],dudz[k,i,j]],[dvdx[k,i,j],dvdy[k,i,j],dvdz[k,i,j]],[dwdx[k,i,j],dwdy[k,i,j],dwdz[k,i,j]]])
            S = 0.5*(dF+dF.T)
            s1[k,i,j] = np.linalg.eig(S)[0].min()
            
    

        