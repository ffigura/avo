# -*- coding: utf-8 -*-
"""
Spyder Editor


Created on Tue Apr 21 17:26:40 2020

@author: Felipe

"""

import numpy as np

def gassmann(vp1,vs1,rho1,phi,k0,k_f1,rho_f1,k_f2,rho_f2):
    """
    Computes elastic properties (Vp, Vs and density) after Gassmann
    fluid substitution.
    
    Avseth et al., Quantitative seismic interpretation, 2006, Page 182.
    
    Parameters
    ----------
    vp1 : float or array
        Initial P-wave velocity - m/s
    vs1 : float or array
        Initial S-wave velocity - m/s       
    rho1 : float or array
        Initial density - g/cm3
    phi : float or array
        Porosity in fraction. 
    k0 : float or array
        Mineral bulk modulus - GPa
    k_f1 : float or array
        Initial fluid bulk modulus - GPa     
    rho_f1 : float or array
        Initial fluid density - g/cm3
    k_f2 : float or array
        Final fluid bulk modulus - GPa
    rho_f2 : float or array
        Final fluid density - g/cm3

    Returns
    -------
    vp2 : float or array
        Final P-wave velocity - m/s
    vs2 : float or array
        Final S-wave velocity - m/s.        
    rho2 : float or array
        Final density - g/cm3.
    """
    vp1=vp1/1000 #km/s
    vs1=vs1/1000 #km/s
    # step1 
    k_sat1=rho1*(vp1**2-(4/3)*vs1**2)
    mu1=rho1*vs1**2
    #step2 - equation 12 for kdry (Smith et al., 2003) 
    #and eq 1 for from Avseth et. al (2006) - the new k_sat2 
    kdry= (k_sat1*((phi*k0)/k_f1+1-phi)-k0)/ \
          ((phi*k0)/k_f1+(k_sat1/k0)-1-phi)
    k_sat2 = kdry+(1- (kdry/k0))**2 /((phi/k_f2) \
             + ((1-phi)/k0) - (kdry/k0**2))
    #step 3 - mu1=mu2
    mu2=mu1
    #step4 - new density
    rho2=rho1+phi*(rho_f2-rho_f1)
    #step 5 - new velocities
    vp2=np.sqrt((k_sat2+(4/3)*mu2)/rho2)*1000 #m/s
    vs2=np.sqrt(mu2/rho2)*1000 #m/s
    
    return(vp2,vs2,rho2)

def vrh(volumes,k,mu):
    #modified from https://github.com/seg/tutorials-2015/blob/master/1506_Seismic_petrophysics_2/Seismic_petrophysics_2.ipynb
    f = np.array(volumes).T
    k = np.resize(np.array(k),np.shape(f))
    mu = np.resize(np.array(mu),np.shape(f))

    k_u = np.sum(f*k, axis=1) #Voigt bound
    k_l = 1. / np.sum(f/k, axis=1) #Reuss bound
    
    mu_u = np.sum(f*mu, axis=1) #Voigt bound
    mu_l = 1. / np.sum(f/mu, axis=1) #Reuss bound
    
    k0 = (k_u+k_l) / 2.  #Hill average
    mu0 = (mu_u+mu_l) / 2.  #Hill average
    
    return(k_u, k_l, mu_u, mu_l, k0, mu0)

def shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    """
    Computes the P-wave reflectivity with Shuey (1985)
    2 and 3 terms for a two-layerd model.
    Avseth et al., Quantitative seismic interpretation,
    2006, Page 182.
    
    Parameters
    ----------
    vp1 : array
        P-wave in the upper layer.
    vs1 : array
        S-wave in the upper layer.
    rho1 : array
        Density in the upper layer.        
    vp2 : array
        P-wave in the lower layer.        
    vs2 : array
        S-wave in the lower layer.      
    rho2 : array
        Density in the lower layer.        
    theta1 : array
        Angles of incidence.

    Returns
    -------
    R0 : array
        Intercept.
    G : array
        Gradient.        
    R2 : array
        Reflection coefficient for the 2-term approximation.
    R3 : array
        Reflection coefficient for the 3-term approximation.   
    """    
    
    theta1 = np.radians(theta1)    
    
    dvp = vp2-vp1
    dvs = vs2-vs1
    drho = rho2-rho1
    vp  = (vp1+vp2)/2
    vs  = (vs1+vs2)/2
    rho = (rho1+rho2)/2   
    
    R0 = 0.5*(dvp/vp + drho/rho)
    G  = 0.5*(dvp/vp) - 2*(vs**2/vp**2)*(drho/rho+2*(dvs/vs))
    F =  0.5*(dvp/vp)
    
    R2 = R0 + G*np.sin(theta1)**2
    R3 = R0 + G*np.sin(theta1)**2 + F*(np.tan(theta1)**2- \
                                       np.sin(theta1)**2)
    
    return (R0,G,R2, R3)

def ai(vp,rho):
    """
    Computes the acoustic impedance

    Parameters
    ----------
    vp : array
        P-velocity.
    rho : array
        Density.

    Returns
    -------
    ai : array
        Acoustic impedance.
    """
    
    ai = vp*rho
    
    return (ai)

def ei(vp,vs,rho,theta1):
    """
    Computes the elasttic impedance.
    Connolly, P., 1999, Elastic impedance: The Leading Edge, 18, 438–452.

    Parameters
    ----------
    vp : array
        P-velocity.
    vs : array
        S-velocity.
    rho : array
        Density.
    theta1 : float
        Incidence angle.

    Returns
    -------
    ei : array
        Elastic impedance.

    """
    theta1 = np.radians(theta1)
    k = (vs/vp)**2
    
    a = 1 + np.tan(theta1)**2
    b = -8*k*np.sin(theta1)**2
    c = 1 - 4*k*np.sin(theta1)**2
    
    ei = (vp**a)*(vs**b)*(rho**c)
    
    
    return (ei)

def pr(vp,vs):
    """
    Computes the Poisson ratio

    Parameters
    ----------
    vp : array
        P-velocity.
    vs : array
        S-velocity.

    Returns
    -------
    pr : array
        Poisson ratio.
    """
    vpvs=vp/vs
    
    pr = 0.5*((vpvs**2-2)/(vpvs**2-1))
    
    return (pr)    