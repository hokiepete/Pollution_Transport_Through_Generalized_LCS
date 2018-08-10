# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 19:34:10 2016

@author: pnolan86

source: https://www.youtube.com/watch?v=mlAuOKD1ff8
"""

from netCDF4 import Dataset
#import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.basemap import Basemap
f = open('NC_MetaData.doc', 'w')
ncfile='subset_wrfout_d01_2011-07-01_00_00_00'
#ncfile="hosiendata.nc"
root = Dataset(ncfile,'r') #read the data

#Query number of dimensions
dims = root.dimensions #dictionary
ndims = len(dims) #number of dimensions
f.write("The # of dimensions = "+str(ndims)+"\n")

# f.write the name and length of each dimension
for key in dims:
    f.write("Dimension["+key+"] = "+str(len(dims[key]))+"\n")
    
gattrs = root.ncattrs() #dictionary
ngattrs = len(gattrs) #number of attributes
f.write("The # of attributes = "+str(ngattrs)+"\n")

#f.write Global Attributes
for key in gattrs:
    f.write("Global_Attribute["+key+"] = "+str(getattr(root,key))+"\n")
    #f.write(

vars = root.variables #dictionary, all variables in dataset
nvars = len(vars) #number of variables in dataset
f.write("The # of variables = "+str(nvars)+"\n")

#For each variable query its dimensions and attributes
for var in vars:
    f.write("---------- variable "+var+" ----------\n")
    f.write("shape = "+str(vars[var].shape)+"\n") #dimensions of variable, tuple
    vdims = vars[var].dimensions #vdims is a tuple
    for vd in vdims:
        f.write("dimension["+vd+"] = "+str(len(dims[vd]))+"\n") #f.write length of variable
    vattrs = vars[var].ncattrs() #dictionary
    f.write("number of attributes = "+str(len(vattrs))+"\n")
    for vat in vattrs:
        f.write("attribute["+vat+"] = "+str(getattr(vars[var],vat))+"\n")
    #now just a slice of data
    #a = vars[var][0:23]
    #f.write(a)
f.close()
root.close()
'''
varshape = vars["vgrid2"].shape
latitude=varshape[0]
longitude=varshape[1]

f = open('Lat&Long.dat', 'w')
f.write("---------- LATITUDE AND LONGITUDE ----------\n")
f.write(str(max(vars["vgrid2"][0,:,0]))+"\n") #longitude max
f.write(str(min(vars["vgrid2"][0,:,0]))+"\n") #longitude min
f.write(str(max(vars["vgrid2"][:,0,1]))+"\n") #latitude max
f.write(str(min(vars["vgrid2"][:,0,1]))+"\n") #latitude min
f.write(str(len(vars["vgrid2"][0,:,0]))+"\n") #longitude length
f.write(str(len(vars["vgrid2"][:,0,1]))+"\n") #latitude length
for lat in range(0,latitude):
    f.write("\n\n")
    for lon in range(0,longitude):
        f.write("["+str(vars["vgrid2"][lat,lon,1])+", "+str(vars["vgrid2"][lat,lon,0])+"]")

f.write("\n\n")
f.write(str(vars["vgrid2"][0,:,0]))
f.write("\n\n")
f.write(str(vars["vgrid2"][0,:,1]))
f.write("\n\n")
f.write(str(vars["vgrid2"][:,0,0]))
f.write("\n\n")
f.write(str(vars["vgrid2"][:,0,1]))
f.close()
print("Done")
'''
"""
+getattr(rootgrp,key))
+str(getattr(rootgrp,key)))
f.write("Dimensions")
f.write(dims)
f.write("Variables")
f.write(vars)
f.write("Attributes")
f.write(attrs)
"""