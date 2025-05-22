import numpy as np
import numexpr as ne
import math
import scipy.integrate as integrate
import h5py

#------------Functions required for estimating mass flux--------------#
def sigmax(M):
    R = getR(M) 
    return np.sqrt(np.log(1. + R*M*M/4))


def getR(M):
    alpha = 2.5
    R = 0.5 * ((3-alpha)/(2-alpha)) * ( 1. - M**(2*(2.-alpha)))/(1. - M**(2.*(3.-alpha)))
    return R


def getPm(x, M):
    sigma_x = sigmax(M)
    yy = (x - sigma_x**2/2.)**2/(2.*sigma_x**2)
    pm = np.exp(-yy)/np.sqrt(2.*math.pi*sigma_x**2)
    return pm


def Uax_fAr(a, gamma, x):
    Ua_x =  np.sqrt((gamma*np.exp(-x) -1.)) * np.sqrt((a-1.)/a) #fixed area
    return Ua_x

def Uax_fSA(a, gamma, x):
    Ua_x = ((gamma*np.exp(-x)-1./a)*(a-1.))**(0.5 )   #fixed solid angle
    return Ua_x

def Uax_int(a, gamma, x):
    Ua_x = (gamma*np.exp(-x)*np.log(a)-(1. -1./a))**(0.5 )   #intermediate
    return Ua_x


def fHI_anly_fAr(a, gamma, x, xi):
    prefac = 2. * xi * np.exp(-x)/np.sqrt(gamma*np.exp(-x)-1.)
    fHI = prefac * np.sqrt((a-1.)/a)
    if(fHI>1.): fHI = 1.
    return fHI
    

def dMassFlux_HI(x, a, gamma, xi, M, p):
    pm = getPm(x, M)
    if(p == 0):
        U_ax = Uax_fAr(a, gamma, x)
        fHI  = fHI_anly_fAr(a, gamma, x, xi)
        
    elif(p==1 ):
        U_ax    =  Uax_int(a, gamma, x) 
        fHI     = fHI_integrand(a, gamma, x, xi, 1) 
    elif(p==2):
        U_ax    =  Uax_fSA(a, gamma, x) 
        fHI     = fHI_integrand(a, gamma, x, xi, 2) 
    if(fHI>1.0):
        fHI = 1.
    return (pm  *  fHI)
    
        
def dMassFlux_tot(x, a, gamma, M, p):
    pm = getPm(x, M)
    if(p == 0):
        U_ax = Uax_fAr(a, gamma, x)
    elif(p==1):
        U_ax    =  Uax_int(a, gamma, x) 
    elif(p==2):
        U_ax    =  Uax_fSA(a, gamma, x) 
    return (pm)


def dMass_HI(x, a, gamma, xi, M, p):
    pm = getPm(x, M)
    if(p == 0):
        U_ax = Uax_fAr(a, gamma, x)
        fHI  = fHI_anly_fAr(a, gamma, x, xi)     
    elif(p==1 ):
        U_ax    =  Uax_int(a, gamma, x) 
        fHI     = fHI_integrand(a, gamma, x, xi, 1) 
    elif(p==2):
        U_ax    =  Uax_fSA(a, gamma, x) 
        fHI     = fHI_integrand(a, gamma, x, xi, 2) 
    if(fHI>1.0):
        fHI = 1.
    return (pm  *  fHI/U_ax )
    
        
def dMass_tot(x, a, gamma, M, p):
    pm = getPm(x, M)
    if(p == 0):
        U_ax = Uax_fAr(a, gamma, x)
    elif(p==1):
        U_ax    =  Uax_int(a, gamma, x) 
    elif(p==2):
        U_ax    =  Uax_fSA(a, gamma, x) 
    return (pm/U_ax)

#---------------------------------------------------------------------#

#-----------Getting fHI, tdyn etc-------------------------------------#


def da_Uax_fSA(a, gamma, x, t0):
    U_ax = Uax_fSA(a, gamma, x)  #fixed solid angle
    dt = 1./(U_ax)
    return dt*t0


def da_Uax_fAr(a, gamma, x, t0):
    U_ax = Uax_fAr(a, gamma, x)  #fixed area
    return U_ax


def da_Uax_int(a, gamma, x, t0):
    U_ax = Uax_int(a, gamma, x)  #intermediate case
    dt = 1./(U_ax)
    return dt*t0


def fHI_fAr(a, gamma, x, t0, G0, Sigma0):
    prefac  = 2. * G0 * t0 * np.exp(-x)/Sigma0/np.sqrt(gamma*np.exp(-x) -1.)
    fHI = prefac * np.sqrt(1. - 1./a)
    if(fHI>1.): fHI = 1.
    return fHI


def fHI_integrand(a, gamma, x, xi, p):
    prefac = xi * np.exp(-x) 
    if(p==1):
        Uax = Uax_int(a, gamma, x) 
    elif(p==2):
        Uax = Uax_fSA(a, gamma, x)
        
    integrand = prefac * (a**(p-2))/Uax
    return integrand



def dfHI_da_fAr(a, gamma, x, t0, G0, Sigma0):
    U_ax    = Uax_fAr(a, gamma, x)  
    dfHI_da =  (G0 * t0 *np.exp(-x) /a/a/U_ax/Sigma0)
    return dfHI_da
