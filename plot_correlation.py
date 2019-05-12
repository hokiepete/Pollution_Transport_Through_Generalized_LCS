# -*- coding: utf-8 -*-
"""
Created on Sat May 11 23:58:35 2019

@author: pnola
"""

import numpy as np
with np.load('correlation.npz') as data:
    c1 = data['c1']
    c2 = data['c2']
    c3 = data['c3']
    z = data['z']
    
    
