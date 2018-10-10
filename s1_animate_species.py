# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 18:45:43 2018

@author: pnola
"""

from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import interp
from scipy import interpolate
import mapping_functions as mf
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker as ticker
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})

titlefont = {'fontsize':12}
labelfont = {'fontsize':10}
tickfont = {'fontsize':8}

grid_spacing = 12*1000 #km 2 m

tdim = 24
xdim = 102
ydim = 82

figwidth = 6
FigSize=(figwidth, ydim/xdim*figwidth)

for species in ['NO2','SO2','O3','ANH4J','ASO4J','ANAJ']:
    root = Dataset('wrf_species_2011_07_01','r')
    cen_lat = getattr(root,'CEN_LAT')
    cen_lon = getattr(root,'CEN_LON')
    true_lat1 = getattr(root,'TRUELAT1')
    true_lat2 = getattr(root,'TRUELAT2')
    ref_lat = getattr(root,'MOAD_CEN_LAT')
    ref_lon = getattr(root,'STAND_LON')
    vars = root.variables
    #Species, Vertically Integrated
    u = vars['U'+species][:,:,:]
    v = vars['V'+species][:,:,:]
    lat = vars['XLAT'][0,:,:]
    lon = vars['XLONG'][0,:,:]
    root.close()
    checklon, checklat = mf.lonlat2km(ref_lon,ref_lat,lon,lat,true_lat1,true_lat2) 
    '''
    if species == '' :
        print('yes')
        u = mf.unstagger(u[:,:,:],2)
        v = mf.unstagger(v[:,:,:],1)
    '''
    if u.shape!=v.shape:
        import sys
        sys.exit('Need to unstagger velocity fields')
    
    x = np.linspace(0,grid_spacing*(xdim-1),xdim)
    dx = x[1]-x[0]
    y = np.linspace(0,grid_spacing*(ydim-1),ydim)
    dy = y[1]-y[0]
    x, y = np.meshgrid(x,y)
    
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
                resolution = 'h',
                area_thresh=1000.,
                )
    parallels = np.arange(round(lat_min,0),lat_max+2,2)
    meridians = np.arange(round(lon_max,0),lon_min-2,-2)
    s1 = np.ma.empty([tdim,ydim,xdim])
    for t in range(u.shape[0]):
        dudy,dudx = np.gradient(u[t,:,:],dy,dx)
        dvdy,dvdx = np.gradient(v[t,:,:],dy,dx)
        J = np.array([[0, 1], [-1, 0]])
        for i in range(ydim):
            for j in range(xdim):
                if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[t,i,j] and v[t,i,j]) is not np.ma.masked:    
                    Utemp = np.array([u[t,i,j], v[t,i,j]])
                    Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
                    S = 0.5*(Grad + np.transpose(Grad))
                    s1[t,i,j] = 3600*np.min(np.linalg.eig(S)[0])
                else:
                    s1[t,i,j] = np.ma.masked
    colormax = 3/4*np.max(-s1,axis=None)
    colormin = 3/4*np.min(-s1,axis=None)
    for t in range(u.shape[0]):
        fig = plt.figure(figsize=FigSize)
        cs = m.contourf(lon,lat,-s1[t,:,:],levels=np.linspace(np.min(-s1[t,:,:],axis=None),np.max(-s1[t,:,:],axis=None),301),latlon=True,vmin=colormin,vmax=colormax)
        m.drawcoastlines()
        m.drawstates()
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
        plt.yticks(**tickfont)
        plt.xticks(**tickfont)
        #plt.title("Time Step = {0:02d} hrs".format(t),fontsize=18)
        plt.title("July 1st 2011, {0:02d}00 GMT".format(t),**titlefont)
        plt.savefig(species+'_S1_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
        plt.close('all')
'''
if species == '' :
            plt.savefig('Wind_S1_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
        else:
            plt.savefig(species+'_S1_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
'''