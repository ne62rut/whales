# -*- coding: utf-8 -*-
"""
Created on Mon May 29 2023

@author: ardhuin
"""

"""waveform_models.py: A Python module for LRM waveform models
"""
#================================================================
# Imports
#----------------------------------------------------------------
import numpy as np
import scipy
from scipy import special

def  myfun_brown_LS(incognita,data)  :
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


def  myfun_brown_ML(incognita,data)  :
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
     ratio = ydata/fff 
     cy= ( ratio - log(ratio)).sum()
     
     return cy


