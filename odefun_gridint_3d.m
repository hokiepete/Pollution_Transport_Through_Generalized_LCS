function [ dydt ] = odefun_gridint_3d( t,Y,U,V,W )
%ODEFUN_GRIDINT_3D Summary of this function goes here
%   Detailed explanation goes here

dydt = zeros(3,1);
dydt(1) = U(Y(1),Y(2),Y(3),t);%interp3(x,y,time,u,Y(1),Y(2),t,'spline');
dydt(2) = V(Y(1),Y(2),Y(3),t);%interp3(x,y,time,v,Y(1),Y(2),t,'spline');
dydt(3) = w(Y(1),Y(2),Y(3),t);%interp3(x,y,time,v,Y(1),Y(2),t,'spline');

end

