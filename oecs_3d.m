clear all
close all
clc
load oecs_data

for i=1:5
    dirdiv1 = smooth3(dirdiv1);
    concav1 = smooth3(concav1);
    s1 = smooth3(s1);
end
dirdiv1(s1>0)=NaN;
dirdiv1(concav1<=0)=NaN;
%dirdiv1(concav1<=0)=NaN;
figure
vol3d_v2('cdata',dirdiv1(:,:,4:end))
colorbar
figure
isosurface(x(:,:,4:end),y(:,:,4:end),z(:,:,4:end),dirdiv1(:,:,4:end),0,'verbose')
axis tight
camlight 
lighting gouraud