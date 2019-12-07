import netCDF4 as nc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
from mapping_functions import lonlat2m
from scipy.interpolate import RegularGridInterpolator
os.environ["PROJ_LIB"] = "C:\\ProgramData\\Anaconda3\\Library\\share"
from mpl_toolkits.basemap import Basemap
from scipy.integrate import solve_ivp
from numpy import nan, isnan
from scipy.interpolate import interp1d

root = nc.Dataset('redux_species_data_Q'+'.nc','r')
vars = root.variables
lon = vars['lon'][:]
lat = vars['lat'][:]
root.close()
lon_min = np.min(lon,axis=None)
lon_max = np.max(lon,axis=None)
lat_min = np.min(lat,axis=None)
lat_max = np.max(lat,axis=None)

m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'l',
            area_thresh=1000.,
            suppress_ticks=True,
            )

xdim = 405
ydim = 325

spec = ['Q','SO2','O3','ASO4J','NO2','ANAJ']

s = spec[5]

root = nc.Dataset('redux_species_data_'+s+'.nc','r')
vars = root.variables
if s=='Q':
    u = vars['eastward_vel'][:,:,:]
    v = vars['northward_vel'][:,:,:]

else:
    u = vars['eastward_vel'][:,:,:]
    v = vars['northward_vel'][:,:,:]


lon = vars['lon'][:]
lat = vars['lat'][:]

t = np.linspace(0,23,24)*3600


"""
tt, tlon, tlat = np.meshgrid(t,lon,lat)
tt = np.transpose(tt,[1,2,0])
tlon = np.transpose(tlon,[1,2,0])
tlat = np.transpose(tlat,[1,2,0])
"""

lon,lat = np.meshgrid(lon,lat)
xx,yy = m(lon,lat)
points = (t,yy[:,0],xx[0,:])

u[u.data==999.]=np.nan
v[v.data==999.]=np.nan

global fu
global fv

fu = RegularGridInterpolator(points,u,fill_value=np.nan)
fv = RegularGridInterpolator(points,v,fill_value=np.nan)

time = 15
with np.load(f'{s}_{time}_att_rep.npz') as file:
    s1_directional_derivative=file['arr_0']
    s2_directional_derivative=file['arr_1']
    s1_concavity=file['arr_2']
    s2_concavity=file['arr_3']
    s1=file['arr_4']
    s2=file['arr_5']
percen = 85

s1_directional_derivative = np.ma.masked_where(s1>=0,s1_directional_derivative)
s2_directional_derivative = np.ma.masked_where(s2<=0,s2_directional_derivative)
s1_directional_derivative = np.ma.masked_where(s1>=-np.percentile(-s1,percen),s1_directional_derivative)
s2_directional_derivative = np.ma.masked_where(s2<=np.percentile(s2,percen),s2_directional_derivative)

#Now we apply the third condition, i.e. concavity 
s1_directional_derivative = np.ma.masked_where(s1_concavity<=0,s1_directional_derivative)
s2_directional_derivative = np.ma.masked_where(s2_concavity>=0,s2_directional_derivative)

#m.contour(lon,lat,s2_directional_derivative,latlon=True,levels=[0],colors='red',linewidths=1.0)
#m.contour(lon,lat,s1_directional_derivative,latlon=True,levels=[0],colors='blue',linewidths=1.0)

test = os.listdir()

for item in test:
    if item.endswith(".png"):
        os.remove(os.path.join(item))
aridge = m.contour(xx,yy,s1_directional_derivative,levels=0,colors='blue')
rridge = m.contour(xx,yy,s2_directional_derivative,levels=0,colors='red')
m.drawcoastlines()
plt.savefig('{:03d}_all.png'.format(000))
plt.close('all')

def fun(t,y):
    Y = y[0]
    X = y[1]
    #print(t,X,Y)
    #input('print')
    if isnan(Y) or isnan(X):
        #return -999., -999. 
        return np.nan, np.nan
    elif y[0]<=0 or y[1]<=0:
        
        return np.nan, np.nan
    else:
        #print(fu((t,Y,X)), fv((t,Y,X)))
        #print(t,Y,X)
        return fu((t,Y,X)), fv((t,Y,X))
    

def event_nan(t, y): 
    #print('y = ',y)
    #if  y[0] == -999. or y[1] == -999.:
    if isnan(y[0]) or isnan(y[1]):
        return 0
    else:
        return 1
event_nan.terminal = True
event_nan.direction = 0

def event_zero(t, y):
    if y[0]<=0 or y[1]<=0:
        return 0
    else:
        return 1

event_nan.terminal = True
event_nan.direction = 0
event_zero.terminal = True
event_zero.direction = 0

tspan = [3600*15,3600*21]
teval = [3600*15,3600*18,3600*21]

pp = aridge.collections[1].get_paths()
flow_y_15 = []
flow_y_18 = []
flow_y_21 = []
flow_x_15 = []
flow_x_18 = []
flow_x_21 = []

for p in pp:#range(len(pp)):
    v = p.vertices#pp[p].vertices
    xx = v[:,0]
    yy = v[:,1]
    line_y_15 = []
    line_y_18 = []
    line_y_21 = []
    line_x_15 = []
    line_x_18 = []
    line_x_21 = []
    for y0 in zip(yy,xx):
        dydt = solve_ivp(fun=fun,t_span=tspan,y0=y0,t_eval=teval,events=[event_nan,event_zero])
        #fy = interp1d(dydt.t, dydt.y[0],fill_value=nan)
        #fx = interp1d(dydt.t, dydt.y[1],fill_value=nan)
        #yt = fy(teval)
        #xt = fx(teval)
        yt = list(dydt.y[0])
        while len(yt) < 3:
            yt.append(nan)
        xt = list(dydt.y[1])
        while len(xt) < 3:
            xt.append(nan)
            
        line_y_15.append(yt[0])
        line_y_18.append(yt[1])
        line_y_21.append(yt[2])
        line_x_15.append(xt[0])
        line_x_18.append(xt[1])
        line_x_21.append(xt[2])
    flow_y_15.append(line_y_15)
    flow_y_18.append(line_y_18)
    flow_y_21.append(line_y_21)
    flow_x_15.append(line_x_15)
    flow_x_18.append(line_x_18)
    flow_x_21.append(line_x_21)
    
            
attracting_ridges = [
        flow_y_15,
        flow_y_18,
        flow_y_21,
        flow_x_15,
        flow_x_18,
        flow_x_21]
    



pp = rridge.collections[1].get_paths()
flow_y_15 = []
flow_y_18 = []
flow_y_21 = []
flow_x_15 = []
flow_x_18 = []
flow_x_21 = []

for p in pp:#range(len(pp)):
    v = p.vertices#pp[p].vertices
    xx = v[:,0]
    yy = v[:,1]
    line_y_15 = []
    line_y_18 = []
    line_y_21 = []
    line_x_15 = []
    line_x_18 = []
    line_x_21 = []
    for y0 in zip(yy,xx):
        dydt = solve_ivp(fun=fun,t_span=tspan,y0=y0,t_eval=teval,events=[event_nan,event_zero])
        #fy = interp1d(dydt.t, dydt.y[0],fill_value=nan)
        #fx = interp1d(dydt.t, dydt.y[1],fill_value=nan)
        #yt = fy(teval)
        #xt = fx(teval)
        yt = list(dydt.y[0])
        while len(yt) < 3:
            yt.append(nan)
        xt = list(dydt.y[1])
        while len(xt) < 3:
            xt.append(nan)
            
        line_y_15.append(yt[0])
        line_y_18.append(yt[1])
        line_y_21.append(yt[2])
        line_x_15.append(xt[0])
        line_x_18.append(xt[1])
        line_x_21.append(xt[2])
    flow_y_15.append(line_y_15)
    flow_y_18.append(line_y_18)
    flow_y_21.append(line_y_21)
    flow_x_15.append(line_x_15)
    flow_x_18.append(line_x_18)
    flow_x_21.append(line_x_21)
    
            
repelling_ridges = [
        flow_y_15,
        flow_y_18,
        flow_y_21,
        flow_x_15,
        flow_x_18,
        flow_x_21]

np.savez(f'ridges_flow_{s}.npz',attracting_ridges=attracting_ridges,repelling_ridges=repelling_ridges)

