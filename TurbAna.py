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
    Purpose: class of mean flow field
        
    Parameters
    ----------
    MeanFlow: 2D numpy array, float.
              nrow = ngrid (number of grid)
              ncol = nvar (e.g., p, T, u, v, w, mut)           
    '''
    
    def __init__(self, MeanFlow=np.ones([1,6])):
        try:
            [ngrid, nvar] = MeanFlow.shape
            self.ngrid = ngrid
            self.nvar  = nvar
            
            if nvar == 4: # 2D data with p, T, u, v
                print('Initialize MeanFlowField using 2D data...')
                MeanFlow = np.concatenate((MeanFlow, np.zeros(shape=[ngrid, 2])), axis=1)
            elif nvar == 5: # 2D data with p, T, u, v, mut
                print('Initialize MeanFlowField using 2D data...')
                MeanFlow = np.concatenate((MeanFlow, np.zeros(shape=[ngrid, 1])), axis=1)                
            elif nvar == 6: # 3D data with p, T, u, v, w, mut
                print('Initialize MeanFlowField using 3D data...')
            else:
                raise ValueError
                
        except:
            raise ValueError('Wrong data type: please use 2D numpy array with ncol = 4 or 9 for MeanFlowField')

        # initialize values
        self.p = MeanFlow[:,0]
        self.T = MeanFlow[:,1] 
        self.Velocity = MeanFlow[:,2:5]
        self.EddyVisc = MeanFlow[:,5]
        self.LamVisc  = self.CalcLamVisc
        self.Density  = self.CalcDensity
        
    def __str__(self):
        return("Mean flow field of %.i grid points."%(self.ngrid))
    
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
        '''
        Purpose: calculate air density using the ideal gas equation
        '''
        R=287        
        rho = self.p/R/self.T
        return(rho)


class MeanGradField:
    '''
    Purpose: class of mean gradient field
        
    Parameters
    ----------
    MeanGrad: 2D numpy array, float.
              nrow = ngrid (number of grid)
              ncol = nvar (e.g., dudx, dudy, dvdx, dvdy for 2D;
                           dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz for 3D)
    '''
    
    def __init__(self, MeanGrad=np.ones([1,9])):
        try:
            [ngrid, nvar] = MeanGrad.shape
            self.ngrid = ngrid
            self.nvar  = nvar
            
            if nvar == 4: # 2D data with uu,vv,ww,uv
                print('Initialize MeanGradField using 2D data...')
                MeanGrad = np.concatenate((MeanGrad[:,0:2],np.zeros(shape=[ngrid, 1]),
                                           MeanGrad[:,2:4],np.zeros(shape=[ngrid, 4])), axis=1)
            elif nvar == 9: # 3D data with 
                print('Initialize MeanGradField using 3D data...')
            else:
                raise ValueError
                
        except:
            raise ValueError('Wrong data type: please use 2D numpy array with ncol = 4 or 9 for MeanGradField')

        # initialize values
        self.VelGrad = MeanGrad
        self.StrainRate = self.CalcStrainRate
        self.Vorticity  = self.CalcVorticity
        
    def __str__(self):
        return("Mean gradient field of %.i grid points."%(self.ngrid))

    @property
    def CalcStrainRate(self):
        '''
        Purpose: calculate strain rate tensor components (Sxx, Syy, Szz, Sxy, Sxz, Syz)
        '''
        Sxx = self.VelGrad[:,0]
        Syy = self.VelGrad[:,4]
        Szz = self.VelGrad[:,8]
        Sxy = 0.5*(self.VelGrad[:,1]+self.VelGrad[:,3])
        Sxz = 0.5*(self.VelGrad[:,2]+self.VelGrad[:,6])
        Syz = 0.5*(self.VelGrad[:,5]+self.VelGrad[:,7])
        S   = (np.array([Sxx,Syy,Szz,Sxy,Sxz,Syz])).T
        return(S)

    @property
    def CalcVorticity(self):
        '''
        Purpose: calculate vorticity tensor components (Oxy, Oxz, Oyz)
        '''
        Oxy = 0.5*(self.VelGrad[:,1]-self.VelGrad[:,3])
        Oxz = 0.5*(self.VelGrad[:,2]-self.VelGrad[:,6])
        Oyz = 0.5*(self.VelGrad[:,5]-self.VelGrad[:,7])
        O   = (np.array([Oxy,Oxz,Oyz])).T
        return(O)
    
class ReynoldsStressTensor:
    '''
    Purpose: class of Reynolds stress tensor
        
    Parameters
    ----------
    tau_ij: 2D numpy array, float.
            nrow = ngrid (number of grid)
            ncol = nvar (e.g., uu,vv,ww,uv for 2D; uu,vv,ww,uv,uw,vw for 3D)           
    '''
    
    def __init__(self, tau_ij=np.ones([1,6])):
        try:
            [ngrid, nvar] = tau_ij.shape
            self.ngrid = ngrid
            self.nvar  = nvar
            
            if nvar == 4: # 2D data with uu,vv,ww,uv
                print('Initialize ReynoldsStressTensor using 2D data...')
                tau_ij = np.concatenate((tau_ij, np.zeros(shape=[ngrid, 2])), axis=1)
            elif nvar == 6: # 3D data with 
                print('Initialize ReynoldsStressTensor using 3D data...')
            else:
                raise ValueError
                
        except:
            raise ValueError('Wrong data type: please use 2D numpy array with ncol = 4 or 6 for ReynoldsStressTensor')
            
        self.Components = tau_ij
        self.TKE = self.CalcTKE
        self.AnisotropyEigenVal = self.CalcAnisotropyEigenVal
        
    def __str__(self):
        return("Reynolds stress tensors of %.i grid points."%(self.ngrid))
    
    @property
    def CalcTKE(self):
        '''
        Purpose: calculate turbulence kinetic energy
        '''        
        k = np.sum(self.Components[:,0:3], axis=1)/2
        return(k)
    
    @property
    def CalcAnisotropyTensor(self):
        '''
        Purpose: calculate turbulence anisotropy tensor
        '''
        k = self.TKE
        k[k<=0] = 1e-9 # avoid division by zero
        a_ij = self.Components/2/np.reshape(k, newshape=(-1,1))
        a_ij[:,0:3] -= 1/3
        return(a_ij)
    
    @property
    def CalcAnisotropyEigenVal(self):
        '''
        Purpose: calculate eigenvalues of turbulence anisotropy tensor
        '''
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
        '''
        Purpose: calculate turbulence data coordinates in Lumley triangle
        '''        
        EigenVals = self.AnisotropyEigenVal
        PC2 = EigenVals[:,0]**2 + EigenVals[:,1]**2 + EigenVals[:,0]*EigenVals[:,1]
        PC3 = -EigenVals[:,0]*EigenVals[:,1]*(EigenVals[:,0] + EigenVals[:,1])
        PCs = np.concatenate((PC3.reshape([-1,1]), PC2.reshape([-1,1])), axis=1)
        return(PCs)
        
    def TurbTriCoor(self):
        '''
        Purpose: calculate turbulence data coordinates in turbulence triangle
        '''  
        EigenVals = self.AnisotropyEigenVal
        xi  = cubic_root((-EigenVals[:,0]*EigenVals[:,1]*(EigenVals[:,0]+EigenVals[:,1]))/2)
        eta = ((EigenVals[:,0]**2 + EigenVals[:,1]**2 + EigenVals[:,0]*EigenVals[:,1])/3)**0.5
        corrs = np.concatenate((xi.reshape([-1,1]), eta.reshape([-1,1])), axis=1)
        return(corrs)
        
    def BaryTriCoor(self):
        '''
        Purpose: calculate turbulence data coordinates in barycentric map
        '''  
        EigenVals = self.AnisotropyEigenVal
        xB = EigenVals[:,0]-EigenVals[:,1]+(3*EigenVals[:,2]+1)/2
        yB = (3*EigenVals[:,2]+1)*np.sqrt(3)/2
        coors = np.concatenate((xB.reshape([-1,1]), yB.reshape([-1,1])), axis=1)
        return(coors)

    def AniRGB(self,c_off=0.65,c_exp=5):
        '''
        Purpose: calculate RGB values from eigenvalues of turbulence anisotropy tensor
        '''  
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

def calc_EddyVisc(RST, MeanGrad, S_ref=100):
    '''
    Purpose: Calculate eddy viscosity based on Boussinesq assumption using 
             Reynolds stress components and velocity gradients
    
    Parameters
    ----------
    RST: ReynoldsStressTensor Class.
    MeanGrad: MeanGradField Class.
    
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
    dudx = MeanGrad.VelGrad[:,0]
    dudy = MeanGrad.VelGrad[:,1]
    dudz = MeanGrad.VelGrad[:,2]
    dvdx = MeanGrad.VelGrad[:,3]
    dvdy = MeanGrad.VelGrad[:,4]
    dvdz = MeanGrad.VelGrad[:,5]
    dwdx = MeanGrad.VelGrad[:,6]
    dwdy = MeanGrad.VelGrad[:,7]
    dwdz = MeanGrad.VelGrad[:,8]
    flag = np.zeros(shape=(ng,1))
    
    # calc eddy viscosity
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


def plot_Lumley_tri(method='w/o data', coors=None):
    '''
    Purpose: plot Lumley triangle template
    
    Parameters
    ----------
    method : string; w/o data (default) for plotting frame only; else for plotting data
    coors  : 2D numpy array, float; coordinates calculated from ReynoldsStressTensor.LumleyTriCoor
             nrow = ngrid (number of grid)
             ncol = 2; Lumley triangle coordinates PC3 and PC2
        
    Returns
    -------
    fig    : matplotlib figure
    '''

    fig = plt.figure(figsize=(6,6))
    
    # plot data
    if method != 'w/o data':
        try:
            plt.scatter(coors[:,0], coors[:,1], zorder=0)
        except:
            raise ValueError("User input 'coors' is not 2D numpy array")
    
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

def plot_turb_tri(method='w/o data', coors=None):
    '''
    Purpose: plot turbulence triangle
    
    Parameters
    ----------
    method : string; w/o data (default) for plotting frame only; else for plotting data
    coors  : 2D numpy array, float; coordinates calculated from ReynoldsStressTensor.TurbTriCoor
             nrow = ngrid (number of grid)
             ncol = 2; turbulence triangle coordinates xi and eta
        
    Returns
    -------
    fig    : matplotlib figure
    '''

    fig = plt.figure(figsize=(6,6))

    # plot data
    if method != 'w/o data':
        try:
            plt.scatter(coors[:,0], coors[:,1], zorder=0)
        except:
            raise ValueError("User input 'coors' is not 2D numpy array")
            
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

def plot_bary_tri(method='w/o data', coors=None):
    '''
    Purpose: plot barycentric map

    Parameters
    ----------
    method : string; w/o data (default) for plotting frame only; else for plotting data
    coors  : 2D numpy array, float; coordinates calculated from ReynoldsStressTensor.BaryTriCoor
             nrow = ngrid (number of grid)
             ncol = 2; barycentric map coordinates xB and yB
        
    Returns
    -------
    fig    : matplotlib figure
    '''

    fig = plt.figure(figsize=(7,6))

    # plot data
    if method != 'w/o data':
        try:
            plt.scatter(coors[:,0], coors[:,1], zorder=0)
        except:
            raise ValueError("User input 'coors' is not 2D numpy array")
    
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

    Parameters
    ----------
    [c_off, c_exp] : float; parameters define the colormap; recommended [0.65, 5] or [0.80, 5]
    nsample        : number of scatter points uniformly distributed in the triangle; default 60
        
    Returns
    -------
    fig    : matplotlib figure

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