# -*- coding: utf-8 -*-
"""
Spyder Editor


Created on Tue Apr 21 17:26:40 2020

@author: Felipe

Least squares regularization with the L1 and L2 norms
"""

import numpy as np

def l2_norm(A,y):
    """
    Least squares or L2-norm solution.
    Aster, R.C., Borchers, B., Thurber, C.H., 2018.
    Parameter estimation and inverse problems, 3rd edition, Page 26

    Parameters
    ----------
    A : 2D array
        The sensitivity matrix.
    y : array
        Input data.

    Returns
    -------
    pest : array
        Estimated parameters
    predict : array
       Predicted data.
    res : array
       Residuals of the fit.
    """

    pest = np.linalg.solve(np.dot(A.T,A),np.dot(A.T,y))
    predict = np.dot(A,pest)
    res = y - predict

    return(pest,predict,res)

def l1_norm(A,y,itmax=20):
    """
    L1-norm solution - Iteratively reweighted least squares IRLS.
    Aster, R.C., Borchers, B., Thurber, C.H., 2018.
    Parameter estimation and inverse problems, 3rd edition, Page 46

    Parameters
    ----------
    A : 2D array
        The sensitivity matrix.
    y : array
        Input data.
    itmax : integer
        The number of iterations.

    Returns
    -------
    pest0 : array
        Estimated parameters.
    predict : array
       Predicted data.
    r : array
       Residuals of the fit. 
    """
    
    #initial guess from L2 norm
    pest0 = np.linalg.solve(np.dot(A.T,A),np.dot(A.T,y))
    predict0 = np.dot(A,pest0)
    r = y - predict0
    phi0 = np.sum(np.abs(r))
    
    for i in range(itmax):
        
        #the process is updated because of the vector r
        R = np.diag(1./np.abs(r+1e-15))
        ATR = np.dot(A.T,R)
        pest = np.linalg.solve(np.dot(ATR,A),np.dot(ATR,y))
        predict = np.dot(A,pest)
        r = y - predict        
        phi = np.sum(np.abs(r))                
        #stop criteria     
        if (np.abs(phi - phi0)/np.abs(phi)) > 1e-5:
            pest0 = np.copy(pest)
            phi0 = np.copy(phi)
            predict0 = np.copy(predict)
        else:
            break              

    return pest0, predict, r
