close all
clear all
clc
load wrfdata_3d
q=1
warning('off','all')
xwant = 405
ywant= 325
zwant = length(z)
%[t,z,y,x] = ndgrid(time,z,y,x);
yy=linspace(0,972000,ywant);
xx=linspace(0,1212000,xwant);
[xx,yy,zz] = meshgrid(xx,yy,z);
[x,y,z,t] = ndgrid(x,y,z,time);
a = [4,3,2,1];
u = permute(u,a);
v = permute(v,a);
w = permute(w,a);
U = griddedInterpolant(x,y,z,t,u,'spline','none')
V = griddedInterpolant(x,y,z,t,v,'spline','none')
W = griddedInterpolant(x,y,z,t,w,'spline','none')

t0 = 21*3600;
tf = 15*3600;
tdim = 37;
t_want = linspace(t0,tf,tdim);

fx = NaN(ywant,xwant,zwant,tdim);
fy = NaN(ywant,xwant,zwant,tdim);
fz = NaN(ywant,xwant,zwant,tdim);

TSPAN = [t0 tf]; % Solve from t=1 to t=5
opts=odeset('event',@eventfun_gridint_wrf_3d,'RelTol',1e-6)%,'AbsTol',1e-14);
for k=(7*(q-1)+1):(7*(q-1)+7)
    sprintf('%03d, %03d, %03d',i,j,k)
    tic
    for i=1:ywant
        for j=1:xwant
            %tic
            Y0=[xx(i,j,k),yy(i,j,k),zz(i,j,k)];
            %[T Y] = ode45(@odefun_gridint, TSPAN, Y0, opts, U,V); % Solve ODE
            [T Y] = ode45(@odefun_gridint_3d, TSPAN, Y0, opts, U,V,W); % Solve ODE
            if sum(~isnan(Y(:,1)))<2
                fx(i,j,k,1) = Y(1,1);
            else
                fx(i,j,k,:) = interp1(T,squeeze(Y(:,1)),t_want,'spline',nan);
            end
            if sum(~isnan(Y(:,2)))<2
                fy(i,j,k,1) = Y(1,2);
            else
                fy(i,j,k,:) = interp1(T,squeeze(Y(:,2)),t_want,'spline',nan);
            end
            if sum(~isnan(Y(:,3)))<2
                fz(i,j,k,1) = Y(1,3);
            else
                fz(i,j,k,:) = interp1(T,squeeze(Y(:,3)),t_want,'spline',nan);
            end


            %fx(i,j,k,:) = interp1(T,squeeze(Y(:,1)),t_want,'spline',nan);
            %fy(i,j,k,:) = interp1(T,squeeze(Y(:,2)),t_want,'spline',nan);
            %fz(i,j,k,:) = interp1(T,squeeze(Y(:,3)),t_want,'spline',nan);
            %toc
        end
    end
    toc
end
time = t_want;
save flow_map_gridint_3d_q1 fx fy fz xx yy zz time

%{

function dydt = odefun_gridint(t,Y,U,V)
    dydt = zeros(2,1);
    dydt(1) = U(Y(1),Y(2),t);%interp3(x,y,time,u,Y(1),Y(2),t,'spline');
    dydt(2) = V(Y(1),Y(2),t);%interp3(x,y,time,v,Y(1),Y(2),t,'spline');
end

function [value,isterminal,direction]=eventfun_gridint(t,Y,U,V)
    value=[Y(1),Y(1)-1212000,Y(2),Y(2)-972000];
    isterminal=[1,1,1,1];
    direction=[0,0,0,0];
end
%}
