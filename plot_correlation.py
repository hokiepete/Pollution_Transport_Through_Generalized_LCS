# -*- coding: utf-8 -*-
"""
Created on Sat May 11 23:58:35 2019

@author: pnola
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
plt.rc('font', **{'family': 'serif', 'serif': ['cmr10']})
titlefont = {'fontsize':12}
labelfont = {'fontsize':12}
tickfont = {'fontsize':10}

with np.load('correlation.npz') as data:
    c1 = data['c1']
    c2 = data['c2']
    c3 = data['c3']
    z = data['z']
    
x=np.linspace(0,34,35)
height = 4.3
siz = [height,3.9/5*height]
list_ = ['h2o','so2','anaj','no2','o3','aso4j','anh4j']
for h in range(len(c1)):
    plt.close('all')
    plt.figure(figsize=siz)
    #plt.subplot(131)
    plt.plot(c1[h],z,label='T=0 hr')
    #plt.subplot(132)
    plt.plot(c2[h],z,label='T=-1 hr')
    #plt.subplot(133)
    plt.plot(c3[h],z,label='T=-6 hr')
    plt.legend(**labelfont, fancybox=True, framealpha=0)
    plt.xlabel('Correlation',**labelfont)
    plt.ylabel('Height AGL (m)',**labelfont)
    plt.xticks(**tickfont)
    plt.yticks(**tickfont)
    plt.savefig('corr_'+list_[h]+'.eps',transparent=False, bbox_inches='tight',pad_inches=0.03)
    
"""
for h in range(len(c1)):
    plt.close('all')
    plt.figure()
    #plt.subplot(131)
    plt.plot(c1[h][:31],z[:31],label='T=0 hr')
    #plt.subplot(132)
    plt.plot(c2[h][:31],z[:31],label='T=-1 hr')
    #plt.subplot(133)
    plt.plot(c3[h][:31],z[:31],label='T=-6 hr')
    plt.legend()
    plt.savefig('corr'+list_[h]+'.png')
    
"""