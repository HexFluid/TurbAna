# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 10:50:40 2021

@author: hxi29
"""

import h5py
import os
import numpy as np

current_path = os.getcwd() # assuming Python launched in the 'TurbAna' dir
data_path    = os.path.join(current_path,'tutorials','bstep_data','bstep_DNS.h5')

h5f  = h5py.File(data_path,'r')
Grid = h5f['grid'][:]         # grid point coordinates
TurbStat = h5f['TurbStat'][:] # Reynolds stress components
h5f.close()

sel_idx = (Grid[:,0]==4)&(Grid[:,1]<=1) # selected profile at x/H=4

# fill in zero for 3D
TurbStat = np.concatenate((TurbStat, np.zeros(shape=[TurbStat.shape[0], 2])), axis=1)



import TurbAna
RST = TurbAna.ReynoldsStressTensor(TurbStat[sel_idx,:])


# plot template
import matplotlib.pyplot as plt
coors = RST.BaryTriCoor()
RGB = RST.AniRGB()
fig = TurbAna.plot_bary_tri()
plt.scatter(coors[:,0], coors[:,1], facecolors=RGB, zorder=0)