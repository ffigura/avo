# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:54:03 2020

@author: Felipe
"""
import numpy as np

def ai(vp,rho):
    """
    Computes de acoustic impedance

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

def nei(vp,vs,rho,vp0,vs0,rho0,theta1):
    """
    Whitcombe, D, 2002, Elastic impedance normalization, Geophysics, 67 (1),
    60–62.

    Parameters
    ----------
    vp : array
        P-velocity.
    vs : array
        S-velocity.
    rho : array
        Density.
    vp0 : float
        Shale reference P-velocity.
    vs0 : float
        Shale reference S-velocity.
    rho0 : float
        Shale reference density.
    theta1 : float
        Incidence angle.

    Returns
    -------
    nei : array
        Normalized elastic impedance.

    """
    
    theta1 = np.radians(theta1)
    k = (vs/vp)**2
    
    a = 1 + np.tan(theta1)**2
    b = -8*k*np.sin(theta1)**2
    c = 1 - 4*k*np.sin(theta1)**2
    
    nei = vp0*rho0*((vp/vp0)**a)*((vs/vs0)**b)*((rho/rho0)**c)
    
    return (nei)