clear all
close all
clc
ncfile='wrf_2011_07_01'
dx = 12000;
dy = 12000;
dimx = 102;
dimy = 82;
dimz = 35
u = permute(unstagger(double(ncread(ncfile,'U')),0),[2,1,3,4]);
v = permute(unstagger(double(ncread(ncfile,'V')),1),[2,1,3,4]);
w = permute(unstagger(double(ncread(ncfile,'W')),2),[2,1,3,4]);
grd = permute(double(ncread(ncfile,'HGT')),[2,1,3,4]);
hgt_in = permute(unstagger(double((ncread(ncfile,'PH') + ncread(ncfile,'PHB')))/9.81,2),[2,1,3,4]);
for h=1:dimz
    hgt(:,:,h)=hgt_in(:,:,h,1)-grd(:,:,1);
end

%look at hgt < 500
%clear hgt_in grd
x = linspace(0,dx*(dimx-1),dimx);
y = linspace(0,dy*(dimy-1),dimy);
z = repmat(0,1,dimz);
[xx,yy,zz] = meshgrid(x,y,z);
dz = 100
z = [1000:dz:4000];
dimz = length(z)
[x,y,z] = meshgrid(x,y,z);
u=u(:,:,:,1);
v=v(:,:,:,1);
w=w(:,:,:,1);
%
u = griddata(xx,yy,hgt,u,x,y,z,'linear');
v = griddata(xx,yy,hgt,v,x,y,z,'linear');
w = griddata(xx,yy,hgt,w,x,y,z,'linear');

t=1;
[dudx,dudy,dudz] = gradient(u(:,:,:),dx,dy,dz);
[dvdx,dvdy,dvdz] = gradient(v(:,:,:),dx,dy,dz);
[dwdx,dwdy,dwdz] = gradient(w(:,:,:),dx,dy,dz);
for i = 1:dimy
    for j = 1:dimx
        for k = 1:dimz
            dF = [dudx(i,j,k),dudy(i,j,k),dudz(i,j,k);dvdx(i,j,k),dvdy(i,j,k),dvdz(i,j,k);dwdx(i,j,k),dwdy(i,j,k),dwdz(i,j,k)];
            S = 0.5*(dF+dF');
            [V,D] = eig(S);
            if ~issorted(diag(D))
                [D,I] = sort(diag(D));
                V = V(:, I);
            end
            s1(i,j,k) = D(1,1);
            s2(i,j,k) = D(end,end);
            X1(i,j,k,:) = V(:,1);
            X2(i,j,k,:) = V(:,end);
            %[V,D] = eig(S);
            %[D,I] = sort(diag(D));
            %V = V(:, I);
            %s1(i,j,k) = min(D);
            %eig1(i,j,k,:) = V(:,1);
        end
    end
end
%
[ds1dx,ds1dy,ds1dz] = gradient(s1,dx,dy,dz);
[ds1dxdx,ds1dxdy,ds1dxdz] = gradient(ds1dx,dx,dy,dz);
[ds1dydx,ds1dydy,ds1dydz] = gradient(ds1dy,dx,dy,dz);
[ds1dzdx,ds1dzdy,ds1dzdz] = gradient(ds1dz,dx,dy,dz);

[ds2dx,ds2dy,ds2dz] = gradient(s2,dx,dy,dz);
[ds2dxdx,ds2dxdy,ds2dxdz] = gradient(ds2dx,dx,dy,dz);
[ds2dydx,ds2dydy,ds2dydz] = gradient(ds2dy,dx,dy,dz);
[ds2dzdx,ds2dzdy,ds2dzdz] = gradient(ds2dz,dx,dy,dz);

%

for i = 1:dimy
    for j=1:dimx
        for k =1:dimz
            concav1(i,j,k) = dot([ds1dxdx(i,j,k),ds1dxdy(i,j,k),ds1dxdz(i,j,k);ds1dydx(i,j,k),ds1dydy(i,j,k),ds1dydz(i,j,k);ds1dzdx(i,j,k),ds1dzdy(i,j,k),ds1dzdz(i,j,k)]*squeeze(X1(i,j,k,:)),squeeze(X1(i,j,k,:)));
            dirdiv1(i,j,k) = dot([ds1dx(i,j,k),ds1dy(i,j,k),ds1dz(i,j,k)],squeeze(X1(i,j,k,:)));
            concav2(i,j,k) = dot([ds2dxdx(i,j,k),ds2dxdy(i,j,k),ds2dxdz(i,j,k);ds2dydx(i,j,k),ds2dydy(i,j,k),ds2dydz(i,j,k);ds2dzdx(i,j,k),ds2dzdy(i,j,k),ds2dzdz(i,j,k)]*squeeze(X2(i,j,k,:)),squeeze(X2(i,j,k,:)));
            dirdiv2(i,j,k) = dot([ds2dx(i,j,k),ds2dy(i,j,k),ds2dz(i,j,k)],squeeze(X2(i,j,k,:)));
            %{
            if or(concav(i,j,k) < 0,s1>-0.5)
                dirdiv(i,j,k) = dot([dsdx(i,j,k),dsdy(i,j,k),dsdz(i,j,k)],squeeze(X1(i,j,k,:)));
            else
                dirdiv(i,j,k) = nan;
            end
            %}
        end
    end
end
%
%dirdiv = dsdx.*X1(:,:,:,1) + dsdy.*X1(:,:,:,2) + dsdz.*X1(:,:,:,3);
%[y,x,z] = meshgrid(y,x,hgt_mean);
save oecs_data s1 s2 dirdiv1 dirdiv2 concav1 concav2 x y z
%figure
%{
isosurface(y,x,z,dirdiv,0)
'isosurface'
FV=isosurface(x,y,z,dirdiv,0,'verbose')
save rawisosurface FV x y z


%}

%{
figure
plot.cdata = dirdiv;
plot.xdata = [0, max(max(max(y)))];
plot.ydata = [0, max(max(max(x)))];
plot.zdata = [hgt_mean(1),hgt_mean(end)];
plot.parent = [];
plot.alpha = [];
plot.texture = '3D'
vol3d_v2(plot)

%}



%z = linspace(0,34,35);

%[y,x,z] = meshgrid(y,x,z);


%}