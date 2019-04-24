close all
clear all
clc

load wrfdata.mat

[t,z,y,x] = ndgrid(time,z,y,x);

U = griddedInterpolant(t,z,y,x,u,'spline','none')
V = griddedInterpolant(t,z,y,x,v,'spline','none')