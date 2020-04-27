# -*- coding: utf-8 -*-
"""
Spyder Editor


Created on Tue Apr 21 17:26:40 2020

@author: Felipe

Approximations of Aki and Richard and Shuey
"""

import numpy as np

def snell(vp1, vp2, vs1, vs2, theta1):
    """
    Calculates the angles of and refraction and reflection for an incident
    P-wave in a two-layered system.
    AVO - Chopra and Castagna, 2014, Page 6.

    :1 is for upper layer and 2 is for lower layer
    :param vp: Compressional velocity
    :param vs: Shear velocity
    :param theta1: Angle of incidence of P-wave
    """
    p = np.sin(theta1)/vp1        # Ray parameter
    thetas1 = np.arcsin(p*vs1)    # S-wave reflection

    # P refraction below first critical angle
    theta2 = np.arcsin(p*vp2)

    # S refraction below second critical angle
    thetas2 = np.arcsin(p*vs2)

    return(theta2, thetas1, thetas2, p)

def akirichards(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    '''
    Calculates P-wave reflectivity with Aki and Richard's (1980) equation
    AVO - Chopra and Castagna, 2014, Page 62.
    
    :1 is for upper layer and 2 is for lower layer
    :param vp: Compressional velocity
    :param vs: Shear velocity
    :param rho: Density
    :param theta1: Angle of incidence of P-wave
    '''
    
    theta1 = np.radians(theta1)    
    theta2, _, _, p = snell(vp1, vp2, vs1, vs2, theta1)
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
    
    return (R1+R2-R3)

def shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    '''
    Calculates P-wave reflectivity with Shuey (1985) 2 and 3 terms
    Avseth et al., Quantitative seismic interpretation, 2006, Page 182.
    
    :1 is for upper layer and 2 is for lower layer
    :param vp: Compressional velocity
    :param vs: Shear velocity
    :param rho: Density
    :param theta1: Angle of incidence of P-wave
    '''    
    
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

def shueyrc(vp0, vs0, rho0, theta):
    '''
    Shuey (1985)
    Calculates P-wave reflectivity with Shuey 2 and 3 terms
    Avseth et al., Quantitative seismic interpretation, 2006, Page 182.
    This implementation gets the whole log
    
    :param vp: Compressional velocity
    :param vs: Shear velocity
    :param rho: Density
    :param theta1: Angle of incidence of P-wave
    '''    
    
    theta = np.radians(theta)
    
    drho=rho0[1:]-rho0[:-1]
    drho=np.insert(drho,0,drho[0]) #insert in the first position
    
    dvp=vp0[1:]-vp0[:-1]
    dvp=np.insert(dvp,0,dvp[0])    
    
    dvs=vs0[1:]-vs0[:-1]
    dvs=np.insert(dvs,0,dvs[0])     
    
    rho=(rho0[1:]+rho0[:-1])/2.0
    rho=np.insert(rho,0,rho[0])

    vp=(vp0[1:]+vp0[:-1])/2.0
    vp=np.insert(vp,0,vp[0])
    
    vs=(vs0[1:]+vs0[:-1])/2.0
    vs=np.insert(vs,0,vs[0])    

    # Compute three-term reflectivity
    r0 = 0.5 * (dvp/vp + drho/rho)
    g = 0.5 * dvp/vp - 2 * (vs**2/vp**2) * (drho/rho + 2 * dvs/vs)
    f = 0.5 * dvp/vp

    term1 = np.outer(r0,1)
    term2 = np.outer(g, np.sin(theta)**2)
    term3 = np.outer(f, (np.tan(theta)**2 - np.sin(theta)**2))
    
    reflection = term1 + term2 + term3
    return (reflection,r0,g)

def rickerwave(f = 25, length = 0.512, dt = 0.004):
    #Ricker wavelet
    time = np.arange(-length/2, (length-dt)/2, dt)
    amplitude = (1.0 - 2.0*(np.pi**2)*(f**2)*(time**2))* \
                          np.exp(-(np.pi**2)*(f**2)*(time**2))
    return time, amplitude

def reflect_coef(ip):
    #Reflection coefficient for plane vertical wave
    rc=(ip[1:]-ip[:-1])/(ip[1:]+ip[:-1])
    rc=np.append(rc,rc[-1])
    
    return(rc)
