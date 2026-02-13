# -*- coding: utf-8 -*-
"""
Created on Mon May 29 2023

@author: M. Passaro & F. ardhuin
"""

"""waveform_models.py: A Python module for LRM waveform models
"""
#================================================================
# Imports
#----------------------------------------------------------------
import numpy as np
import scipy
from scipy import special

######################  Defines waveform theoretical models: brown from WHALES code, includes PTR convolution
def  wf_brown_eval(xdata,incognita,noise,Gamma,Zeta,c_xi,PTR)  :
    '''
    Define a waveform as in the WHALES code.
    Based on a Brown model.
    inputs :
            - xdata : range gates
            - incognita : (3,) vector with [0] = epoch, [1] = Hs, [2] = amplitude
            - noise : thermal noise
            - Gamma : coeff with antenna bandwidth (gamma in Tourain et al. 2020)
            - Zeta : off nadir pointing angle (= mispointing)
            - c_xi : 4 c /(G * h) (with G = Gamma in Tourain et al. 2020)
            - PTR : optionnal PTR
    output : - waveform
    '''

    ff0=noise+( incognita[2]/2*np.exp((-4/Gamma)*(np.sin(Zeta))**2) \
    * np.exp (-  c_xi*( (xdata-incognita[0])-c_xi*incognita[1]**2/2) ) \
    *   (  1+scipy.special.erf( ((xdata-incognita[0])-c_xi*incognita[1]**2)/((np.sqrt(2)*incognita[1]))  ) ) \
    )
    if PTR[0] < 1:
       fff =fftconvolve(ff0,PTR,mode='same')
    else:
       fff=ff0

    return fff
    
######################  Cost functions for retracking ... should use function above. 
def  waveform_brown_LS(incognita,data)  :
     """
     returns the least-square distance between the waveform data[0] and the theoretical 
     Brown-Hayne functional form, The unknown parameters in this version (17 Dec 2013) are Epoch, Sigma and Amplitude, where 
     sigma=( sqrt( (incognita(2)/(2*0.3)) ^2+SigmaP^2) ) is the rising time of the leading edge
     
     For the explanation of the terms in the equation, please check "Coastal Altimetry" Book
     
     """
     
     ydata =data[0] #Waveform coefficients
     Gamma =data[1]
     Zeta  =data[2]
     xdata =data[3]  #Epoch
     SigmaP=data[4]
     c_xi  =data[5]  #Term related to the slope of the trailing edge
     weights=data[6]  #Weights to apply to the residuals
         
     fff = ( incognita[2]/2*np.exp((-4/Gamma)*(np.sin(Zeta))**2) \
     * np.exp (-  c_xi*( (xdata-incognita[0])-c_xi*incognita[1]**2/2) ) \
     *   (  1+scipy.special.erf( ((xdata-incognita[0])-c_xi*incognita[1]**2)/((np.sqrt(2)*incognita[1]))  ) ) \
     )
    
     cy= (   weights *  ((ydata - fff) **2)).sum()
     
     return cy


def  waveform_brown_ML(incognita,data)  :
     """
     returns the ML distance between the waveform data[0] and the theoretical 
     Brown-Hayne functional form, The unknown parameters in this version (17 Dec 2013) are Epoch, Sigma and Amplitude, where 
     sigma=( sqrt( (incognita(2)/(2*0.3)) ^2+SigmaP^2) ) is the rising time of the leading edge
     
     For the explanation of the terms in the equation, please check "Coastal Altimetry" Book
     
     """
     
     ydata =data[0] #Waveform coefficients
     Gamma =data[1]
     Zeta  =data[2]
     xdata =data[3]  #Epoch
     SigmaP=data[4]
     c_xi  =data[5]  #Term related to the slope of the trailing edge
     weights=data[6]  #Weights to apply to the residuals
         
     fff = ( incognita[2]/2*np.exp((-4/Gamma)*(np.sin(Zeta))**2) \
     * np.exp (-  c_xi*( (xdata-incognita[0])-c_xi*incognita[1]**2/2) ) \
     *   (  1+scipy.special.erf( ((xdata-incognita[0])-c_xi*incognita[1]**2)/((np.sqrt(2)*incognita[1]))  ) ) \
     )
     ratio = np.divide(ydata+1.e-5,fff+1.e-5) 
     cy= ( ratio - np.log(ratio)).sum()
     
     return cy

def brown_model(x, data):
    """
    Compute the theoretical Brown-Hayne model waveform for a given parameter
    vector x = [epoch, sigma, amplitude].

    This function extracts all physical/geophysical parameters from args_tuple
    and evaluates the Brown functional form (leading edge + trailing edge)
    used in classical ocean altimetry retracking.

    Parameters
    ----------
    x : array_like, shape (3,)
        Model parameters:
        - x[0]: epoch (in nanoseconds)
        - x[1]: sigma (leading edge rise time)
        - x[2]: amplitude of the waveform

    data : tuple
        Contains the waveform data and auxiliary constants:
        (ydata, Gamma, Zeta, xdata, SigmaP, c_xi, weights, weightflag)

    Returns
    -------
    fff : ndarray
        The modeled Brown-Hayne waveform evaluated at xdata.
    """    
    ydata = data[0]
    Gamma = data[1]
    Zeta  = data[2]
    xdata = data[3]
    SigmaP= data[4]
    c_xi  = data[5]

    # Build the Brown model waveform fff
    fff = ( x[2]/2*np.exp((-4/Gamma)*(np.sin(Zeta))**2)
        * np.exp(-c_xi * ((xdata - x[0]) - c_xi * x[1]**2 / 2))
        * (1 + scipy.special.erf(((xdata - x[0]) - c_xi * x[1]**2)
                                 / (np.sqrt(2)*x[1])))
    )
    return fff

def brown_residuals(x, data, weights, weightflag):
    """
    Residual vector for Levenberg-Marquardt or Gauss-Newton optimization.

    Computes the pointwise difference between the observed waveform ydata
    and the theoretical Brown-Hayne model waveform produced by brown_model(),
    optionally applying user-supplied weights.

    This function returns the residual vector required by
    scipy.optimize.least_squares().

    Parameters
    ----------
    x : array_like, shape (3,)
        Current estimate of the model parameters.

    data : tuple
        Same tuple passed to brown_model():
        (ydata, Gamma, Zeta, xdata, SigmaP, c_xi, weights, weightflag)

    Returns
    -------
    residuals : ndarray
        The weighted or unweighted residuals: (ydata - model).
    """
#    ydata = data[0]
#    fff   = brown_model(x, data)
#    resid = ydata - fff
#    return resid * weights if weightflag else resid

    ydata = data[0]

    # -------------------------------------------
    # SAFETY 1 — sigma must be positive & non-zero
    # -------------------------------------------
    if x[1] <= 0:
        return np.full_like(ydata, 1e10)

    # -------------------------------------------
    # Compute model waveform
    # -------------------------------------------
    fff = brown_model(x, data)

    # -------------------------------------------
    # SAFETY 2 — model must be finite
    # -------------------------------------------
    if not np.all(np.isfinite(fff)):
        return np.full_like(ydata, 1e10)

    # -------------------------------------------
    # Residual computation
    # -------------------------------------------
    resid = ydata - fff

    return np.sqrt(weights) * resid




