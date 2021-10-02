"""
Turbulence Analyzer (TurbAna) Python toolkit Ver 1.0

The script was originally written for analyzing the turbulence anisotropy of 
compressor tip leakage flows:
    
  He, X., Zhao, F., & Vahdati, M. (2022). Detached Eddy Simulation: Recent
  Development and Application to Compressor Tip Leakage Flow. ASME Journal
  of Turbomachinery, 144(1), 011009.
  
An explict reference to the work above is highly appreciated if this script is
useful for your research.  

Xiao He (xiao.he2014@imperial.ac.uk)
Last update: 30-Sep-2021
"""

# -------------------------------------------------------------------------
# import libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

# -------------------------------------------------------------------------
# main classes
class MeanFlowField:
    '''
    Parameters
    ----------
    MeanFlow: 2D numpy array, float.
              nrow = ngrid (number of grid)
              ncol = nvar (e.g., p, T, u, v, w, mut)
    MeanGrad: 2D numpy array, float.
              nrow = ngrid (number of grid)
              ncol = nvar (e.g., dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)              
    '''
    def __init__(self, MeanFlow=np.ones([1,6]), MeanGrad=np.ones([1,9])):
        self.p = MeanFlow[:,0]
        self.T = MeanFlow[:,1] 
        self.Velocity = MeanFlow[:,2:5]
        self.EddyVisc = MeanFlow[:,5]
        self.LamVisc  = self.CalcLamVisc
        self.Density  = self.CalcDensity
        self.MeanGrad = MeanGrad
        
    def __str__(self):
        return("Mean flow field of %.i grid points."%(self.MeanFlow.shape[0]))
    
    @property
    def CalcLamVisc(self):
        '''
        Purpose: calculate laminar viscosity using Sutherland's law
        '''
        mu_ref = 1.716 * 10**(-5)
        Tc =273.15
        S = 110.4
        mu = mu_ref*(self.T/Tc)**1.5*(Tc+S)/(self.T+S)
        return(mu)
    
    @property
    def CalcDensity(self):
        R=287        
        rho = self.p/R/self.T
        return(rho)
    
    
class ReynoldsStressTensor:
    '''
    Parameters
    ----------
    tau_ij: 2D numpy array, float.
            nrow = ngrid (number of grid)
            ncol = nvar (e.g., uu,vv,ww,uv,uw,vw)           
    '''
    def __init__(self, tau_ij=np.ones([1,6])):
        self.Components = tau_ij
        self.TKE = self.CalcTKE
        self.AnisotropyEigenVal = self.CalcAnisotropyEigenVal
        
    def __str__(self):
        return("Reynolds stress tensors of %.i grid points."%(self.Components.shape[0]))
    
    @property
    def CalcTKE(self):
        k = np.sum(self.Components[:,0:3], axis=1)/2
        return(k)
    
    @property
    def CalcAnisotropyTensor(self):
        k = self.TKE
        k[k<=0] = 1e-9 # avoid division by zero
        a_ij = self.Components/2/np.reshape(k, newshape=(-1,1))
        a_ij[:,0:3] -= 1/3
        return(a_ij)
    
    @property
    def CalcAnisotropyEigenVal(self):
        a_ij = self.CalcAnisotropyTensor
        EigenVals = []
        EigenVal_flags = [] # flags used for debug
        
        # loop over each grid point
        for igrid in range(self.Components.shape[0]):
            matrix = np.array([[a_ij[igrid,0],a_ij[igrid,3],a_ij[igrid,4]],
                               [a_ij[igrid,3],a_ij[igrid,1],a_ij[igrid,5]],
                               [a_ij[igrid,4],a_ij[igrid,5],a_ij[igrid,2]]])
            lambdas = linalg.eigvals(matrix) # solve for eigenvalues, i.e., lambdas
            lambdas_real = np.real(lambdas)
            lambdas_real[::-1].sort() # descending order of eigenvalues
            EigenVals.append(lambdas_real)
            
            # check if eigenvalues are correct
            if np.sum(np.imag(lambdas)) > 0.001:
                print('Non-real eigenvalues detected! Please check data.')
                print('Error occurs at index: '+str(igrid))
                EigenVal_flags.append(-1)
            if np.max(np.real(lambdas))>2/3 or np.min(np.real(lambdas))<-2/3:
                EigenVal_flags.append(0)
                print('Eigenvalues beyond range [-2/3,2/3] detected! Please check data.')
                print('Error occurs at index: '+str(igrid))            
            else:
                EigenVal_flags.append(1)
                
        return(np.array(EigenVals))

    def LumleyTriCoor(self):
        EigenVals = self.AnisotropyEigenVal
        PC2 = EigenVals[:,0]**2 + EigenVals[:,1]**2 + EigenVals[:,0]*EigenVals[:,1]
        PC3 = -EigenVals[:,0]*EigenVals[:,1]*(EigenVals[:,0] + EigenVals[:,1])
        PCs = np.concatenate((PC3.reshape([-1,1]), PC2.reshape([-1,1])), axis=1)
        return(PCs)
        
    def TurbTriCoor(self):
        EigenVals = self.AnisotropyEigenVal
        xi  = cubic_root((-EigenVals[:,0]*EigenVals[:,1]*(EigenVals[:,0]+EigenVals[:,1]))/2)
        eta = ((EigenVals[:,0]**2 + EigenVals[:,1]**2 + EigenVals[:,0]*EigenVals[:,1])/3)**0.5
        corrs = np.concatenate((xi.reshape([-1,1]), eta.reshape([-1,1])), axis=1)
        return(corrs)
        
    def BaryTriCoor(self):
        EigenVals = self.AnisotropyEigenVal
        xB = EigenVals[:,0]-EigenVals[:,1]+(3*EigenVals[:,2]+1)/2
        yB = (3*EigenVals[:,2]+1)*np.sqrt(3)/2
        coors = np.concatenate((xB.reshape([-1,1]), yB.reshape([-1,1])), axis=1)
        return(coors)

    def AniRGB(self,c_off=0.65,c_exp=5):
        EigenVals = self.AnisotropyEigenVal
        R = (EigenVals[:,0]-EigenVals[:,1]+c_off)**c_exp
        G = (2*(EigenVals[:,1]-EigenVals[:,2])+c_off)**c_exp
        B = (3*EigenVals[:,2]+1+c_off)**c_exp
        RGB = np.concatenate((R.reshape([-1,1]), G.reshape([-1,1]), B.reshape([-1,1])), axis=1)
        RGB[RGB>1] = 1
        return(RGB)

# -------------------------------------------------------------------------
# sub-functions
def cubic_root(data):
    '''
    Purpose: Calculate cubic root and return the real root
    
    Parameters
    ----------
    data: numpy array, float.
    
    Return
    root: numpy array, float.
    '''
    data_shape=data.shape
    data=data.flatten()
    
    root=[]
    for idx in range(data.shape[0]):
        if data[idx] >= 0:
            root.append(data[idx]**(1/3))
        else:
            root.append(-(-data[idx])**(1/3))
    root=np.array(root)    
    root=root.reshape(data_shape)
    return root

def calc_EddyVisc(RST, MeanFlow, S_ref=100):
    '''
    Purpose: Calculate eddy viscosity based on Boussinesq assumption using 
             Reynolds stress components and velocity gradients
    
    Parameters
    ----------
    RST: ReynoldsStressTensor Class.
    MeanFlow: MeanFlowField Class.
    
    Return
    EddyVisc: 1D numpy array, float. Calculated dynamic eddy viscosity
    flag:     indicator of regions affected by limiters.
              0: no limiter; 1: S_ref applied; -1: EddyVisc=max{0,self} applied.
    '''
    
    # initialize parameters
    ng = (RST.Components).shape[0]
    uu = RST.Components[:,0]
    vv = RST.Components[:,1]
    ww = RST.Components[:,2]
    uv = RST.Components[:,3]
    uw = RST.Components[:,4]
    vw = RST.Components[:,5]
    dudx = MeanFlow.MeanGrad[:,0]
    dudy = MeanFlow.MeanGrad[:,1]
    dudz = MeanFlow.MeanGrad[:,2]
    dvdx = MeanFlow.MeanGrad[:,3]
    dvdy = MeanFlow.MeanGrad[:,4]
    dvdz = MeanFlow.MeanGrad[:,5]
    dwdx = MeanFlow.MeanGrad[:,6]
    dwdy = MeanFlow.MeanGrad[:,7]
    dwdz = MeanFlow.MeanGrad[:,8]
    flag = np.zeros(shape=(ng,1))
    
    denominator = 2*(dudx**2+dvdy**2+dwdz**2+
                     2*(0.5*(dudy+dvdx))**2+
                     2*(0.5*(dudz+dwdx))**2+
                     2*(0.5*(dvdz+dwdy))**2-
                     1/3*(dudx+dvdy+dwdz)**2)
    numerator = -(uu*dudx+vv*dvdy+ww*dwdz+
                  2*uv*0.5*(dudy+dvdx)+
                  2*uw*0.5*(dudz+dwdx)+
                  2*vw*0.5*(dvdz+dwdy)-
                  1/3*(uu+vv+ww)*(dudx+dvdy+dwdz))
    
    # limit |denominator|>= S_ref**2
    S_ref2 = S_ref**2
    idx_sref_neg = (denominator>=-S_ref2)&(denominator<= 0)
    idx_sref_pos = (denominator<= S_ref2)&(denominator>= 0)
    denominator[idx_sref_neg]=-S_ref2
    denominator[idx_sref_pos]= S_ref2
    flag[idx_sref_neg]=1
    flag[idx_sref_pos]=1    

    # calc eddy viscosity    
    EddyVisc = numerator/denominator
    
    # limit EddyVisc >= 0
    idx_negEddyVisc = (EddyVisc<0)
    EddyVisc[idx_negEddyVisc] = 0
    flag[idx_negEddyVisc] = -1
    
    return EddyVisc, flag


def plot_Lumley_tri():
    '''
    Purpose: plot Lumley triangle template
    '''

    fig = plt.figure(figsize=(6,6))
    
    # figure format
    plt.xlabel('III')
    plt.ylabel('II')
    plt.axis([-0.02,0.08,-0.05,0.4])
    plt.yticks(np.linspace(0,0.4,5))
    
    # triangle bounds
    lam_1s = [2/3, 1/6, 0]
    lam_2s = [-1/3, 1/6, 0]
    for iedge in range(3):
        lam_1=np.linspace(lam_1s[iedge],lam_1s[(iedge+1)%3],100)
        lam_2=np.linspace(lam_2s[iedge],lam_2s[(iedge+1)%3],100)
        PC2_temp=lam_1**2+lam_2**2+lam_1*lam_2
        PC3_temp=-lam_1*lam_2*(lam_1+lam_2)
        plt.plot(PC3_temp,PC2_temp,color='grey')
    
    # texts
    plt.text(-0.0027,-0.025,'3C',fontsize=16)
    plt.text(-0.018,0.075,'2C',fontsize=16)
    plt.text(0.07,0.34,'1C',fontsize=16)
    
    return fig

def plot_turb_tri():
    '''
    Purpose: plot turbulence triangle template
    '''

    fig = plt.figure(figsize=(6,6))
    
    # figure format
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$\eta$')
    plt.axis([-0.2,0.4,-0.03,0.35])
    plt.yticks(np.linspace(0,0.3,4))
    
    # triangle bounds
    lam_1s = [2/3, 1/6, 0]
    lam_2s = [-1/3, 1/6, 0]
    for iedge in range(3):
        lam_1=np.linspace(lam_1s[iedge],lam_1s[(iedge+1)%3],100)
        lam_2=np.linspace(lam_2s[iedge],lam_2s[(iedge+1)%3],100)
        eta_temp=((lam_1**2+lam_2**2+lam_1*lam_2)/3)**(1/2)
        xi_temp=cubic_root(-(lam_1*lam_2*(lam_1+lam_2))/2)
        plt.plot(xi_temp,eta_temp,color='grey')    
    
    # texts
    plt.text(-0.02,-0.018,'3C',fontsize=16)
    plt.text(-0.19,0.18,'2C',fontsize=16)
    plt.text(0.34,0.33,'1C',fontsize=16)
    
    return fig

def plot_bary_tri():
    '''
    Purpose: plot barycentric map template
    '''

    fig = plt.figure(figsize=(7,6))
    
    # figure format
    plt.xlabel('$x_B$')
    plt.ylabel('$y_B$')
    plt.axis([-0.1,1.1,-0.01,0.9])
    plt.axis('off')
    
    # triangle bounds
    plt.plot([1,0],[0,0],color='grey')
    plt.plot([1,1/2],[0,np.sqrt(3)/2],color='grey')
    plt.plot([0,1/2],[0,np.sqrt(3)/2],color='grey')
    
    # texts
    plt.text(0.47,0.91,'3C',fontsize=16)
    plt.text(-0.1,-0.05,'2C',fontsize=16)
    plt.text(1.05,-0.05,'1C',fontsize=16)

    return fig

def plot_bary_tri_colormap(c_off=0.65,c_exp=5,nsample=60):
    '''
    Purpose: plot barycentric color map
    ---
    Recommended [c_off, c_exp] values: [0.65, 5]; [0.80, 5]
    '''
    
    fig = plt.figure(figsize=(3.5,3))
    
    # calculate scatter coordinates
    coors=[]
    for irow in range(1,nsample+1):
        l_temp = (irow-1)/(nsample-1)
        x_temp = 0.5*(1-l_temp)
        y_temp = np.sqrt(3)/2*(1-l_temp)
        for icol in range(1,irow+1):
            if irow != 1:
                x_temp2 = x_temp + l_temp*(icol-1)/(irow-1)
            else:
                x_temp2 = x_temp
            coors.append([x_temp2, y_temp])
    coors=np.array(coors)
    
    # calculate scatter facecolor
    c1c = coors[:,0] - coors[:,1]/np.sqrt(3)
    c3c = coors[:,1]*2/np.sqrt(3)
    c2c = 1-c1c-c3c
    R = (c1c+c_off)**c_exp
    G = (c2c+c_off)**c_exp
    B = (c3c+c_off)**c_exp
    RGB=np.concatenate((R.reshape([-1,1]), G.reshape([-1,1]), B.reshape([-1,1])), axis=1)
    RGB=np.array(RGB,dtype=np.float16)
    RGB[RGB>1]=1
    RGB[RGB<0]=0
    
    # plot
    plt.scatter(coors[:,0], coors[:,1], facecolors=RGB, alpha=0.8, s=40, zorder=0)
 
    # figure format
    plt.xlabel('$x_B$')
    plt.ylabel('$y_B$')
    plt.axis([-0.1,1.1,-0.01,0.9])
    plt.axis('off')
    
    # triangle bounds
    plt.plot([1,0],[0,0],color='grey')
    plt.plot([1,1/2],[0,np.sqrt(3)/2],color='grey')
    plt.plot([0,1/2],[0,np.sqrt(3)/2],color='grey')

    # wall boundary
    plt.fill_between([-0.1,0,0.5,0.55], [-0.01,0,np.sqrt(3)/2,np.sqrt(3)/2], 0.9, facecolor='white')
    plt.fill_between([-0.1,1.1], [0,0], -0.1, facecolor='white')
    plt.fill_between([1.1,1,0.5,0.45], [-0.01,0,np.sqrt(3)/2,np.sqrt(3)/2], 0.9, facecolor='white')
    
    # texts
    plt.text(0.43,0.91,'3C',fontsize=16)
    plt.text(-0.16,-0.05,'2C',fontsize=16)
    plt.text(1.05,-0.05,'1C',fontsize=16)
    
    return fig

# End