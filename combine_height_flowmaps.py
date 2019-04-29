from scipy.io import loadmat, savemat
import numpy as np

#fx fy fz xx yy zz time

fx = np.empty((325,405,35,37))
fy = np.empty((325,405,35,37))
fz = np.empty((325,405,35,37))

for i in range(5):
    fx[:,:,7*i:7*i+7,:] = loadmat('flow_map_gridint_3d_q{}.mat'.format(i+1))['fx'][:,:,7*i:7*i+7,:]
    fy[:,:,7*i:7*i+7,:] = loadmat('flow_map_gridint_3d_q{}.mat'.format(i+1))['fy'][:,:,7*i:7*i+7,:]
    fz[:,:,7*i:7*i+7,:] = loadmat('flow_map_gridint_3d_q{}.mat'.format(i+1))['fz'][:,:,7*i:7*i+7,:]

x=loadmat('flow_map_gridint_3d_q1.mat')['xx']
y=loadmat('flow_map_gridint_3d_q1.mat')['yy']
z=loadmat('flow_map_gridint_3d_q1.mat')['zz']
t=loadmat('flow_map_gridint_3d_q1.mat')['time']

savemat('flow_map_wrf_3d.mat',{'fx':fx,'fy':fy,'fz':fz,'x':x,'y':y,'z':z,'t':t})
