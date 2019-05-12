# -*- coding: utf-8 -*-
"""
Created on Sat May 11 23:58:35 2019

@author: pnola
"""
import matplotlib.pyplot as plt
import numpy as np
with np.load('correlation.npz') as data:
    c1 = data['c1']
    c2 = data['c2']
    c3 = data['c3']
    z = data['z']
    
x=np.linspace(0,34,35)
    
for h in range(len(c1)):
    plt.close('all')
    plt.figure()
    plt.subplot(131)
    plt.plot(c1[h],z)
    plt.subplot(132)
    plt.plot(c2[h],z)
    plt.subplot(133)
    plt.plot(c3[h],z)
    plt.savefig('{:02d}.png'.format(h))