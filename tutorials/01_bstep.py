"""
Application of TurbAna to backward-facing step spanwise-averaged and time-averaged
Reynolds stress components data. Details of the data can be found in the following:

  Le, H., Moin, P., & Kim, J. (1997). Direct Numerical Simulation of Turbulent 
  Flow Over a Backward-Facing Step. Journal of Fluid Mechanics, 330, 349-374.
    
  He, X., Zhao, F., & Vahdati, M. (2022). Detached Eddy Simulation: Recent
  Development and Application to Compressor Tip Leakage Flow. ASME Journal
  of Turbomachinery, 144(1), 011009.

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
#                 ncol = 2 (e.g., x/H, y/H)
# TurbStat shape: nrow = ngrid (number of grid)
#                 ncol = nvar (e.g., uu,vv,ww,uv)
# MeanFlow shape: nrow = ngrid (number of grid)
#                 ncol = nvar (e.g., p,T,u,v,mu_t)
# MeanGrad shape: nrow = ngrid (number of grid)
#                 ncol = nvar (e.g., dudx,dudy,dvdx,dvdy)
# -------------------------------------------------------------------------

# Sec. 1 start time
start_sec1 = time.time()

# data path
data_path = os.path.join(current_path, 'bstep_data')

# option to save SPOD results
save_fig  = False  # postprocess figs
save_path = data_path

# load data from h5 format
h5f  = h5py.File(os.path.join(data_path,'bstep_DNS.h5'),'r')
DNS_Grid = h5f['grid'][:]         # grid point coordinates
DNS_TurbStat = h5f['TurbStat'][:] # Reynolds stress components
h5f.close()

nskip = 10
h5f  = h5py.File(os.path.join(data_path,'bstep_DDES.h5'),'r')
DDES_Grid = h5f['grid'][::nskip,:]         # grid point coordinates
DDES_TurbStat = h5f['TurbStat'][::nskip,:] # Reynolds stress components
DDES_MeanFlow = h5f['MeanFlow'][::nskip,:] # Mean flow fields
DDES_MeanGrad = h5f['MeanGrad'][::nskip,:] # Velocity gradients
DDES_Flag     = h5f['FlagDES'][::nskip,:]  # DES flags (shielding function, delta/H, FKH)
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
DNS_RST = TurbAna.ReynoldsStressTensor(DNS_TurbStat)
DDES_RST = TurbAna.ReynoldsStressTensor(DDES_TurbStat)

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
    if legend != 'None':
        plt.legend(loc=legend)

def bstep_frame():
    '''
    Purpose: template for backward-facing step 2D contour plot
    '''

    fig = plt.figure(figsize=(6,4))

    # wall boundary
    plt.fill_between([-4,0,0,30], [1,1,0,0], -1, facecolor='whitesmoke')
    plt.plot([-4,0,0,30], [1,1,0,0],color='black',linewidth=1)

    # figure format
    figure_format('x/H','y/H', [-2,10,-0.5,2],'None')

    return fig

# -------------------------------------------------------------------------   
### 3.1 plot Lumley triangle
# select plot location
xloc = 4
ymax = 1
eps  = 0.01
sel_idx_DNS  = (np.abs(DNS_Grid[:,0]-xloc)<eps)&(DNS_Grid[:,1]<=ymax)
sel_idx_DDES = (np.abs(DDES_Grid[:,0]-xloc)<eps)&(DDES_Grid[:,1]<=ymax)

# plot template
fig31 = TurbAna.plot_Lumley_tri()

# plot DNS data
DNS_coors = DNS_RST.LumleyTriCoor()
plt.scatter(DNS_coors[sel_idx_DNS,0], DNS_coors[sel_idx_DNS,1], label='DNS', 
            facecolor= 'white', edgecolors='red', linewidths = 1.0, zorder = 3)

# plot DDES data
DDES_coors = DDES_RST.LumleyTriCoor()
plt.scatter(DDES_coors[sel_idx_DDES,0], DDES_coors[sel_idx_DDES,1], label='DDES', 
            marker='s', facecolor= 'white', edgecolors='orange', linewidths = 1.0, zorder = 3)

# figure format
plt.legend(loc='upper left')

if save_fig:
    plt.savefig(os.path.join(save_path,'x=%.1f'%xloc+'H_Lumley_tri.png'), dpi=300, bbox_inches='tight')
    plt.close()

print('Plot Lumley triangle finished')

# -------------------------------------------------------------------------   
### 3.2 plot turbulence triangle
# select plot location
xloc = 4
ymax = 1
eps  = 0.01
sel_idx_DNS = (np.abs(DNS_Grid[:,0]-xloc)<eps)&(DNS_Grid[:,1]<=ymax)
sel_idx_DDES = (np.abs(DDES_Grid[:,0]-xloc)<eps)&(DDES_Grid[:,1]<=ymax)

# plot template
fig32 = TurbAna.plot_turb_tri()

# plot DNS data
DNS_coors = DNS_RST.TurbTriCoor()
plt.scatter(DNS_coors[sel_idx_DNS,0], DNS_coors[sel_idx_DNS,1], label='DNS', 
            facecolor= 'white', edgecolors='red', linewidths = 1.0, zorder = 3)

# plot DDES data
DDES_coors = DDES_RST.TurbTriCoor()
plt.scatter(DDES_coors[sel_idx_DDES,0], DDES_coors[sel_idx_DDES,1], label='DDES', 
            marker='s', facecolor= 'white', edgecolors='orange', linewidths = 1.0, zorder = 3)

# figure format
plt.legend(loc='lower right')

if save_fig:
    plt.savefig(os.path.join(save_path,'x=%.1f'%xloc+'H_turb_tri.png'), dpi=300, bbox_inches='tight')
    plt.close()

print('Plot turbulence triangle finished')

# -------------------------------------------------------------------------   
### 3.3 plot barycentric map
# select plot location
xloc = 4
ymax = 1
eps  = 0.01
sel_idx_DNS = (np.abs(DNS_Grid[:,0]-xloc)<eps)&(DNS_Grid[:,1]<=ymax)
sel_idx_DDES = (np.abs(DDES_Grid[:,0]-xloc)<eps)&(DDES_Grid[:,1]<=ymax)

# plot template
fig33 = TurbAna.plot_bary_tri()

# plot DNS data
DNS_coors = DNS_RST.BaryTriCoor()
plt.scatter(DNS_coors[sel_idx_DNS,0], DNS_coors[sel_idx_DNS,1], label='DNS', 
            facecolor= 'white', edgecolors='red', linewidths = 1.0, zorder = 3)

# plot DDES data
DDES_coors = DDES_RST.BaryTriCoor()
plt.scatter(DDES_coors[sel_idx_DDES,0], DDES_coors[sel_idx_DDES,1], label='DDES', 
            marker='s', facecolor= 'white', edgecolors='orange', linewidths = 1.0, zorder = 3)

# figure format
plt.legend(loc='lower center', columnspacing=0.7,ncol=2,bbox_to_anchor=(0.5, -0.15))

if save_fig:
    plt.savefig(os.path.join(save_path,'x=%.1f'%xloc+'H_bary_map.png'), dpi=300, bbox_inches='tight')
    plt.close()

print('Plot barycentric map finished')

# -------------------------------------------------------------------------   
### 3.4 plot anisotropy contour
# colormap option
c_off=0.65
c_exp=5

# calculate RGB value
RGB = DDES_RST.AniRGB(c_off=c_off,c_exp=c_exp)

# plot colormap
fig341 = TurbAna.plot_bary_tri_colormap(c_off=c_off,c_exp=c_exp)
if save_fig:
    plt.savefig(os.path.join(save_path,'anisotropy_colormap.png'), dpi=300, bbox_inches='tight')
    plt.close()

# plot anisotropy 2D contour
fig342 = bstep_frame()
plt.scatter(DDES_Grid[:,0], DDES_Grid[:,1], facecolors=RGB, alpha=0.8, s=10, zorder=0)
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


# -------------------------------------------------------------------------
# 4. Calculate turbulent viscosity
# function TurbAna.calc_EddyVisc(ReynoldsStressTensor, MeanFlowField)
# -------------------------------------------------------------------------

# Sec. 4 start time
start_sec4 = time.time()

# main function
DDES_Mean = TurbAna.MeanFlowField(DDES_MeanFlow)
DDES_Grad = TurbAna.MeanGradField(DDES_MeanGrad)
S_ref = 30
method = 'QCR2013V' # recommend Boussinesq/QCR2013V
[DDES_EddyViscCal, DDES_EddyViscFlag] = TurbAna.calc_EddyVisc(DDES_RST,DDES_Grad,S_ref=S_ref,method=method)

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
fig51 = bstep_frame()
cntr = plt.tricontourf(DDES_Grid[:,0], DDES_Grid[:,1], DDES_EddyViscFlag[:,0], 
                       np.linspace(-0.5,3.5,5),cmap=cm.rainbow,extend='both', zorder=0)

# colorbar
plt.colorbar(cntr,ticks=np.linspace(0,3,4),shrink=0.8,extendfrac='auto',\
             orientation='vertical', label=r'$f_{lim}$')

if save_fig:
    plt.savefig(os.path.join(save_path,'sref=%.i'%S_ref+'_ViscRatio_limiter.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()

# -------------------------------------------------------------------------   
### 5.2 plot calculated eddy/laminar viscosity ratio
fig52 = bstep_frame()
DDES_ViscRatio = DDES_EddyViscCal*DDES_Mean.Density/(DDES_Mean.LamVisc*100) # 100 comes from Re_ref in CFD
cntr = plt.tricontourf(DDES_Grid[:,0], DDES_Grid[:,1], np.log10(1e-3+DDES_ViscRatio), 
                       np.linspace(0,3,21),cmap=cm.gist_ncar,extend='both', zorder=0)

# colorbar
cbar = plt.colorbar(cntr,ticks=np.linspace(0,3,4),shrink=0.8,extendfrac='auto',\
             orientation='vertical', label=r'$\nu_t/\nu$')
cbar.ax.set_yticklabels(['$10^0$', '$10^1$', '$10^2$', '$10^3$'])

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