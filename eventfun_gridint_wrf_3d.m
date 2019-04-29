function [ value,isterminal,direction ] = eventfun_gridint_wrf_3d( t,Y,U,V,W )
%EVENTFUN_GRIDINT_WRF_3D Summary of this function goes here
%   Detailed explanation goes here
    value=[isnan(Y(1))-1,isnan(Y(2))-1,isnan(Y(3))-1];
    isterminal=[1,1,1];
    direction=[0,0,0];


end

