# -*- coding: utf-8 -*-
"""
Spyder Editor


Created on Tue Apr 21 17:26:40 2020

@author: Felipe

Approximations of Aki and Richard and Shuey
"""

import numpy as np

def snell(vp1, vp2, theta1):
    """
    Computes the angles of refraction for an incident P-wave in a two-layered
    model. AVO - Chopra and Castagna, 2014, Page 6.

    Parameters
    ----------
    vp1 : array
        P-wave in the upper layer.
    vp2 : array
        P-wave in the lower layer.
    theta1 : array
        Angles of incidence.

    Returns
    -------
    theta2 : array
        Angles of refraction.
    p : array
       Ray parameter.
    """
    p = np.sin(theta1)/vp1  
    theta2 = np.arcsin(p*vp2)

    return(theta2, p)

def akirichards(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    """
    Computes the P-wave reflectivity with Aki and Richard's (1980) equation for
    a two-layered model.
    AVO - Chopra and Castagna, 2014, Page 62.
    
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
    R  : array
        Reflection coefficient.
    """    
    
    theta1 = np.radians(theta1)    
    theta2, p = snell(vp1, vp2, theta1)
    theta = (theta1 + theta2) / 2.

    dvp = vp2-vp1
    dvs = vs2-vs1
    drho = rho2-rho1
    vp  = (vp1+vp2)/2
    vs  = (vs1+vs2)/2
    rho = (rho1+rho2)/2    
    
    R1 = 0.5*(1-4*p**2*vs**2)*drho/rho
    R2 = 0.5/(np.cos(theta)**2)*dvp/vp
    R3 = 4*p**2*vs**2*dvs/vs
    
    R = R1+R2-R3
    
    return (R)

def shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    """
    Computes the P-wave reflectivity with Shuey (1985) 2 and 3 terms for a 
    two-layerd model.
    Avseth et al., Quantitative seismic interpretation, 2006, Page 182.
    
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
    R3 = R0 + G*np.sin(theta1)**2 + F*(np.tan(theta1)**2-np.sin(theta1)**2)
    
    return (R0,G,R2, R3)

def shueyrc(vp0, vs0, rho0, theta1):
    """
    Computes the P-wave reflectivity with Shuey (1985) 2 terms for a 
    log.
    Avseth et al., Quantitative seismic interpretation, 2006, Page 182.
    
    Parameters
    ----------
    vp0 : array
        P-wave.
    vs0 : array
        S-wave.
    rho0 : array
        Density.        
    theta1 : array
        Angles of incidence.

    Returns
    -------
    R : array
        Reflection coefficient for the 2-term approximation.    
    R0 : array
        Intercept.
    G : array
        Gradient.        
    """      
    
    theta1 = np.radians(theta1)
    
    dvp=vp0[1:]-vp0[:-1]
    dvs=vs0[1:]-vs0[:-1]
    drho=rho0[1:]-rho0[:-1]
    #insert in the first position    
    drho=np.insert(drho,0,drho[0]) 
    dvp=np.insert(dvp,0,dvp[0])    
    dvs=np.insert(dvs,0,dvs[0])     

    vp=(vp0[1:]+vp0[:-1])/2.0
    vs=(vs0[1:]+vs0[:-1])/2.0    
    rho=(rho0[1:]+rho0[:-1])/2.0

    vp=np.insert(vp,0,vp[0])
    vs=np.insert(vs,0,vs[0])    
    rho=np.insert(rho,0,rho[0])

    # Compute two-term reflectivity
    R0 = 0.5 * (dvp/vp + drho/rho)
    G = 0.5 * dvp/vp - 2 * (vs**2/vp**2) * (drho/rho + 2 * dvs/vs)

    term1 = np.outer(R0,1)
    term2 = np.outer(G, np.sin(theta1)**2)
    
    R = term1 + term2 
    return (R,R0,G)

def rickerwave(f = 25, length = 0.512, dt = 0.004):
    """
    Computes the ricker wavelet.
    
    Parameters
    ----------
    f : float
       Central frequency.
    length : float
        Length of the wavelet.
    dt : float
        Sample rate.        

    Returns
    -------
    time : array
        Time.    
    Amplitude : array
        Amplitude of the wavelet.
    """          
    time = np.arange(-length/2, (length-dt)/2, dt)
    amplitude = (1.0 - 2.0*(np.pi**2)*(f**2)*(time**2))* \
                          np.exp(-(np.pi**2)*(f**2)*(time**2))
    return (time, amplitude)

def reflect_coef(ip):
    """
    Computes the reflection coefficient for a plane incident P-wave.
    
    Parameters
    ----------
    ip : array
       P-impedance.

    Returns
    -------
    rc : array
        The reflection coefficient
    """ 
    rc=(ip[1:]-ip[:-1])/(ip[1:]+ip[:-1])
    rc=np.append(rc,rc[-1])
    
    return(rc)
