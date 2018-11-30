clear all
close all
clc
ncfile='wrf_2011_07_01'
dx = 12000;
dy = 12000;
dimx = 102;
dimy = 82;
dimz = 35
u = unstagger(double(ncread(ncfile,'U')),0);
v = unstagger(double(ncread(ncfile,'V')),1);
w = unstagger(double(ncread(ncfile,'W')),2);
grd = double(ncread(ncfile,'HGT'));
hgt_in = unstagger(double((ncread(ncfile,'PH') + ncread(ncfile,'PHB')))/9.81,2);
for t=1:dimz
    hgt(:,:,t)=hgt_in(:,:,t,1)-grd(:,:,1);
end
%clear hgt_in grd
x = linspace(0,dx*(dimx-1),dimx);
y = linspace(0,dy*(dimy-1),dimy);
hgt_mean = squeeze(mean(mean(hgt)));
dz = cat(1,hgt_mean(1), diff(hgt_mean));
%dz = cat(1, diff(hgt_mean),hgt_mean(end));
t=1;
[dudy,dudx,dudz] = gradient(u(:,:,:,t),dy,dx,dz);
[dvdy,dvdx,dvdz] = gradient(v(:,:,:,t),dy,dx,dz);
[dwdy,dwdx,dwdz] = gradient(w(:,:,:,t),dy,dx,dz);
for i = 1:dimx
    for j = 1:dimy
        for k = 1:dimz
            dF = [dudx(i,j,k),dudy(i,j,k),dudz(i,j,k);dvdx(i,j,k),dvdy(i,j,k),dvdz(i,j,k);dwdx(i,j,k),dwdy(i,j,k),dwdz(i,j,k)];
            S = 0.5*(dF+dF');
            [V,D] = eig(S);
            if ~issorted(diag(D))
                [D,I] = sort(diag(D));
                V = V(:, I);
            end
            s1(i,j,k) = D(1,1);
            X1(i,j,k,:) = V(:,1);
            %[V,D] = eig(S);
            %[D,I] = sort(diag(D));
            %V = V(:, I);
            %s1(i,j,k) = min(D);
            %eig1(i,j,k,:) = V(:,1);
        end
    end
end
%
[dsdy,dsdx,dsdz] = gradient(s1,dy,dx,dz);
[dsdydy,dsdydx,dsdydz] = gradient(dsdy,dy,dx,dz);
[dsdxdy,dsdxdx,dsdxdz] = gradient(dsdx,dy,dx,dz);
[dsdzdy,dsdzdx,dsdzdz] = gradient(dsdz,dy,dx,dz);
%

for i = 1:dimx
    for j=1:dimy
        for k =1:dimz
            concav(i,j,k) = dot([dsdxdx(i,j,k),dsdxdy(i,j,k),dsdxdz(i,j,k);dsdydx(i,j,k),dsdydy(i,j,k),dsdydz(i,j,k);dsdzdx(i,j,k),dsdzdy(i,j,k),dsdzdz(i,j,k)]*squeeze(X1(i,j,k,:)),squeeze(X1(i,j,k,:)));
            dirdiv(i,j,k) = dot([dsdx(i,j,k),dsdy(i,j,k),dsdz(i,j,k)],squeeze(X1(i,j,k,:)));
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
[y,x,z] = meshgrid(y,x,hgt_mean);
save oecs_data s1 dirdiv concav x y z
%figure
%isosurface(y,x,z,dirdiv,0)
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


