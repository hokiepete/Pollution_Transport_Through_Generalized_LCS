close all
clear all
clc
load flow_map_wrf_3d
x=double(x);
y=double(y);
t=double(t);
dx = x(1,2,1)-x(1,1,1)
dy = y(2,1,1)-y(1,1,1)
dz = z(1,1,:)-z(1,1,1)

[leny,lenx,lenz,tlen]=size(fx)
T = t-t(1)

sigma = zeros([leny,lenx,lenz,tlen]);
for t = 1:tlen
    T(t)
    %need to skip t0 to avoid division by zero
    %set FTLE to zero instead
    if T(t) == 0
        %sigma(:,:,t) = zeros([leny,lenx])
        continue
    else
        [dfxdx,dfxdy,dfxdz] = gradient(squeeze(fx(:,:,:,t)),dx,dy,dz);
        [dfydx,dfydy,dfydz] = gradient(squeeze(fy(:,:,:,t)),dx,dy,dz);
        [dfzdx,dfzdy,dfzdz] = gradient(squeeze(fz(:,:,:,t)),dx,dy,dz);
        
        for k = 1:lenz
            for i = 1:leny
                for j = 1:lenx
                    if ~isnan(dfxdx(i,j,k))&&~isnan(dfxdy(i,j,k))&&~isnan(dfxdz(i,j,k))&&~isnan(dfydx(i,j,k))&&~isnan(dfydy(i,j,k))&&~isnan(dfydz(i,j,k))&&~isnan(dfzdx(i,j,k))&&~isnan(dfzdy(i,j,k))&&~isnan(dfzdz(i,j,k))
                        gradF = [dfxdx(i,j,k),dfxdy(i,j,k),dfxdz(i,j,k);
                                 dfydx(i,j,k),dfydy(i,j,k),dfydz(i,j,k);
                                 dfzdx(i,j,k),dfzdy(i,j,k),dfzdz(i,j,k)];
                        C = gradF'*gradF;
                        lambda = max(eig(C));
                        sigma(i,j,k,t) = 1/(2*abs(T(t)))*log(lambda);
                    else
                        sigma(i,j,k,t) = nan;
                    end
                end
            end
        end
    end
end

save('FTLE_wrf_3d.mat', 'sigma', 'T','x','y','z');
