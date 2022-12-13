"""
Application of TurbAna to transonic bump pitch-averaged and time-averaged
Reynolds stress components data. Details of the data can be found in the following:

  Uzun, A., & Malik, M. R. (2019). Wall-resolved large-eddy simulations of 
  transonic shock-induced flow separation. AIAA Journal, 57(5), 1955-1972.

Xiao He (xiao.he2014@imperial.ac.uk)
Last update: 30-June-2022
"""

# -------------------------------------------------------------------------
# 0. Import libraries
# standard python libraries
import sys
import time
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import numpy as np
import h5py
import pylab

# TurbAna path
current_path = os.getcwd()
parrent_path = os.path.dirname(current_path)
TurbAna_path = parrent_path

# TurbAna library
sys.path.insert(0, TurbAna_path)
import TurbAna

# -------------------------------------------------------------------------
# 1. Load input data for TurbAna
# grid shape    : nrow = ngrid (number of grid)
#                 ncol = 2 (e.g., x/c, y/c)
# TurbStat shape: nrow = ngrid (number of grid)
#                 ncol = nvar (e.g., uu,vv,ww,uv)
# MeanFlow shape: nrow = ngrid (number of grid)
#                 ncol = nvar (e.g., p,T,u,v)
# MeanGrad shape: nrow = ngrid (number of grid)
#                 ncol = nvar (e.g., dudx,dudy,dvdx,dvdy)
# -------------------------------------------------------------------------

# Sec. 1 start time
start_sec1 = time.time()

# data path
data_path = os.path.join(current_path, 'bump_data')

# option to save TurbAna results
save_fig  = False  # postprocess figs
save_path = data_path

# load data from h5 format
h5f  = h5py.File(os.path.join(data_path,'bump_LES.h5'),'r')
LES_Grid = h5f['grid'][:]         # grid point coordinates
LES_TurbStat = h5f['TurbStat'][:] # Reynolds stress components
LES_MeanFlow = h5f['MeanFlow'][:] # Mean flow fields
LES_MeanGrad = h5f['MeanGrad'][:] # Velocity gradients
h5f.close()

# Sec. 1 end time
end_sec1 = time.time()

print('--------------------------------------'    )
print('Input data loaded!'                        )
print('Time lapsed: %.2f s'%(end_sec1-start_sec1) )
print('--------------------------------------'    )


# -------------------------------------------------------------------------
# 2. Calculate turbulence anisotropy
# function TurbAna.ReynoldsStressTensor(TurbStat)
# -------------------------------------------------------------------------

# Sec. 2 start time
start_sec2 = time.time()

# main function
LES_RST = TurbAna.ReynoldsStressTensor(LES_TurbStat)

# Sec. 2 end time
end_sec2 = time.time()

print('--------------------------------------'    )
print('Calculation of anisotropy finished!'     )
print('Time lapsed: %.2f s'%(end_sec2-start_sec2) )
print('--------------------------------------'    )


# -------------------------------------------------------------------------
# 3. Visualize turbulence anisotropy
# Figs: 1. Lumley triangle
#       2. turbulence triangle
#       3. barycentric map
#       4. anisotropy contour
# -------------------------------------------------------------------------

# Sec. 3 start time
start_sec3 = time.time()

# -------------------------------------------------------------------------
### 3.0 pre-defined function
params={
'axes.labelsize': '20',
'xtick.labelsize': '16',
'ytick.labelsize': '16',
'lines.linewidth': 1.5,
'legend.fontsize': '14',
'figure.figsize': '8, 6'    # set figure size
}
pylab.rcParams.update(params)

def figure_format(xtitle, ytitle, zoom, legend):
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.axis(zoom)
    plt.axis('off')
    if legend != 'None':
        plt.legend(loc=legend)

def bump_frame():
    '''
    Purpose: template for backward-facing step 2D contour plot
    '''

    # fig = plt.figure(figsize=(6,4))
    fig = plt.figure(figsize=(4,2.7))

    # wall boundary
    c=1
    h=1.905/20.32
    r=(h**2+c**2/4)/2/h
    xo=0.5
    yo=0.375-(r-h)
    xlocs = np.concatenate((np.linspace(-2,0,3),np.linspace(0,1,30),np.linspace(1,2,3)))
    ylocs = np.concatenate((np.linspace(0.375,0.375,3),
                            yo+np.sqrt(r**2-(np.linspace(0,1,30)-xo)**2),
                            np.linspace(0.375,0.375,3)))
    plt.fill_between(xlocs, ylocs, 0.35, facecolor='whitesmoke')
    plt.plot(xlocs, ylocs,color='black',linewidth=1)

    # figure format
    figure_format('x/c','y/c', [-0.5,1.5,0.35,0.60],'None')

    return fig

# -------------------------------------------------------------------------   
### 3.1 plot Lumley triangle
# select plot location
xloc = 1
ymin = 0.37501
ymax = 0.55
eps  = 0.005
sel_idx_LES = (np.abs(LES_Grid[:,0]-xloc)<eps)&(LES_Grid[:,1]<=ymax)&(LES_Grid[:,1]>=ymin)

# plot template
fig31 = TurbAna.plot_Lumley_tri()

# plot LES data
LES_coors = LES_RST.LumleyTriCoor()
plt.scatter(LES_coors[sel_idx_LES,0], LES_coors[sel_idx_LES,1], label='LES', 
            marker='s', facecolor= 'white', edgecolors='orange', linewidths = 1.0, zorder = 3)

# figure format
plt.legend(loc='upper left')

if save_fig:
    plt.savefig(os.path.join(save_path,'x=%.1f'%xloc+'c_Lumley_tri.png'), dpi=300, bbox_inches='tight')
    plt.close()

print('Plot Lumley triangle finished')

# -------------------------------------------------------------------------   
### 3.2 plot turbulence triangle
# select plot location
xloc = 1
ymin = 0.37501
ymax = 0.55
eps  = 0.005
sel_idx_LES = (np.abs(LES_Grid[:,0]-xloc)<eps)&(LES_Grid[:,1]<=ymax)&(LES_Grid[:,1]>=ymin)

# plot template
fig32 = TurbAna.plot_turb_tri()

# plot LES data
LES_coors = LES_RST.TurbTriCoor()
plt.scatter(LES_coors[sel_idx_LES,0], LES_coors[sel_idx_LES,1], label='LES', 
            marker='s', facecolor= 'white', edgecolors='orange', linewidths = 1.0, zorder = 3)

# figure format
plt.legend(loc='lower right')

if save_fig:
    plt.savefig(os.path.join(save_path,'x=%.1f'%xloc+'c_turb_tri.png'), dpi=300, bbox_inches='tight')
    plt.close()

print('Plot turbulence triangle finished')

# -------------------------------------------------------------------------   
### 3.3 plot barycentric map
# select plot location
xloc = 1
ymin = 0.37501
ymax = 0.55
eps  = 0.005
sel_idx_LES = (np.abs(LES_Grid[:,0]-xloc)<eps)&(LES_Grid[:,1]<=ymax)&(LES_Grid[:,1]>=ymin)

# plot template
fig33 = TurbAna.plot_bary_tri()

# plot LES data
LES_coors = LES_RST.BaryTriCoor()
plt.scatter(LES_coors[sel_idx_LES,0], LES_coors[sel_idx_LES,1], label='LES', 
            marker='s', facecolor= 'white', edgecolors='orange', linewidths = 1.0, zorder = 3)

# figure format
plt.legend(loc='lower center', columnspacing=0.7,ncol=2,bbox_to_anchor=(0.5, -0.15))

if save_fig:
    plt.savefig(os.path.join(save_path,'x=%.1f'%xloc+'c_bary_map.png'), dpi=300, bbox_inches='tight')
    plt.close()

print('Plot barycentric map finished')

# -------------------------------------------------------------------------   
### 3.4 plot anisotropy contour
# colormap option
c_off=0.65
c_exp=5

# calculate RGB value
RGB = LES_RST.AniRGB(c_off=c_off,c_exp=c_exp)

# plot colormap
fig341 = TurbAna.plot_bary_tri_colormap(c_off=c_off,c_exp=c_exp)
if save_fig:
    plt.savefig(os.path.join(save_path,'anisotropy_colormap.png'), dpi=300, bbox_inches='tight')
    plt.close()

# plot anisotropy 2D contour
fig342 = bump_frame()
plt.scatter(LES_Grid[:,0], LES_Grid[:,1], facecolors=RGB, alpha=0.8, s=10, zorder=0)
if save_fig:
    plt.savefig(os.path.join(save_path,'anisotropy_contour.png'), dpi=300, bbox_inches='tight')
    plt.close()

# -------------------------------------------------------------------------  
# Sec. 3 end time
end_sec3 = time.time()

print('--------------------------------------'    )
print('Anisotropy plots finished!'                )
print('Figs saved to the directory:'              )
print( save_path                                  )
print('Time lapsed: %.2f s'%(end_sec3-start_sec3) )
print('--------------------------------------'    )

#%%
# -------------------------------------------------------------------------
# 4. Calculate turbulent viscosity
# function TurbAna.calc_EddyVisc(ReynoldsStressTensor, MeanFlowField)
# -------------------------------------------------------------------------

# Sec. 4 start time
start_sec4 = time.time()

# main function
LES_Mean = TurbAna.MeanFlowField(LES_MeanFlow)
LES_Grad = TurbAna.MeanGradField(LES_MeanGrad)

S_ref  = 1000
method = 'QCR2013V' # recommend Boussinesq/QCR2013V
[LES_EddyViscCal, LES_EddyViscFlag] = TurbAna.calc_EddyVisc(LES_RST,LES_Grad,S_ref=S_ref,method=method)

# Sec. 4 end time
end_sec4 = time.time()

print('--------------------------------------'     )
print('Calculation of viscosity finished!'         )
print('Time lapsed: %.2f s'%(end_sec4-start_sec4)  )
print('--------------------------------------'     )


# -------------------------------------------------------------------------
# 5. Visualize turbulent viscosity
# Figs: 1. limiter of turbulent viscosity calculation contour
#       2. turbulent viscosity contour
# -------------------------------------------------------------------------

# Sec. 5 start time
start_sec5 = time.time()

# -------------------------------------------------------------------------   
### 5.1 plot eddy viscosity calculation flag
fig51 = bump_frame()
cntr = plt.tricontourf(LES_Grid[:,0], LES_Grid[:,1], LES_EddyViscFlag[:,0], 
                       np.linspace(-0.5,3.5,5),cmap=cm.rainbow,extend='both', zorder=0)

# colorbar
plt.colorbar(cntr,ticks=np.linspace(0,3,4),shrink=0.8,extendfrac='auto',\
             orientation='vertical', label=r'$f_{lim}$')

if save_fig:
    plt.savefig(os.path.join(save_path,'sref=%.i'%S_ref+'_'+method+'_ViscRatio_limiter.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()

# -------------------------------------------------------------------------   
### 5.2 plot calculated eddy/laminar viscosity ratio
fig52 = bump_frame()
LES_ViscRatio = LES_EddyViscCal*LES_Mean.Density/LES_Mean.LamVisc
cntr = plt.tricontourf(LES_Grid[:,0], LES_Grid[:,1], np.log10(1e-3+LES_ViscRatio),
                       np.linspace(0,4,21),cmap=cm.gist_ncar,extend='both', zorder=0)

# colorbar
cbar = plt.colorbar(cntr,ticks=np.linspace(0,4,3),shrink=0.8,extendfrac='auto',\
             orientation='vertical',label=r'$\nu_t/\nu$')
cbar.ax.set_yticklabels(['$10^0$', '$10^2$', '$10^4$'])

if save_fig:
    plt.savefig(os.path.join(save_path,'sref=%.i'%S_ref+'_'+method+'_ViscRatio_contour.png'), dpi=300, bbox_inches='tight')
    plt.close()

# -------------------------------------------------------------------------  
# Sec. 5 end time
end_sec5 = time.time()

print('--------------------------------------'    )
print('Viscosity plots finished!'                 )
print('Figs saved to the directory:'              )
print( save_path                                  )
print('Time lapsed: %.2f s'%(end_sec5-start_sec5) )
print('--------------------------------------'    )    
    
# End