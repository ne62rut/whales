# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 15:16:13 2018

@author: passaro
"""

from scipy.io import netcdf
from scipy import interpolate
from scipy import ndimage
import numpy as np
#import netCDF4
import pickle
from datetime import datetime
from sidereal import *
import os.path
from scipy.io import matlab

def compute_instr_corr_SWH_WHALES(SWH_ALES,my_path_instr_corr_SWH,mission,f='empty'):
    """
    
    28.08.2018 : updated for J3 WHALES
    
    19/01/2016 : the function works with Jason-2    
    
    Inputs:
    1. SWH_ALES = SWH in output of ALES retracker
    2. mypath_instr_corr_SWH = path and name of the instr_corr_SWH look up table 
    3. f = If not inserted, the algorithm computes again the interpolating function        
    
    
    Output:
    
    1. SWH_ALES = new SWH instrumental-corrected
    2. f = to be called in next iteration in order to avoid re-computing the interpolating 2d function
    
    NOTE FOR ENVISAT: even if we use the Sigma0 and the SWH from the SGDR in the Abdalla+look-up table 
    model, we won't obtain the ssb correction that the SGDR use. This is because the Wind value from the 
    altimeter is modified before computing the SSB by using external data.

    The user manual states:''The wind speed is computed (in m/s), using a linear interpolation in the 
    input wind table, according to the algorithm proposed by Abdalla [R17]. The algorithm is based on a 
    fit between EnviSat Ku-band Sigma0 and the collocated ECMWF model wind speed. The result was then 
    adjusted based on in-situ wind measurements.''    
    
    
    """    

         
    
    if isinstance(f, str) == False:    
        mat = matlab.loadmat(my_path_instr_corr_SWH)
        corrvalues = mat['SWHinstrcorr_J3WHALESsgdrD_corrvalues'].squeeze()
        SWHvalues = mat['SWHinstrcorr_J3WHALESsgdrD_swhvalues'].squeeze()           
        
        if SWH_ALES < np.min(SWHvalues):         
            instr_corr_SWH = corrvalues[0]
        elif SWH_ALES > np.max(SWHvalues):           
            instr_corr_SWH = corrvalues[-1]            
        else:
            instr_corr_SWH = f(SWH_ALES)   
                
    else:
        mat = matlab.loadmat(my_path_instr_corr_SWH)
        corrvalues = mat['SWHinstrcorr_J3WHALESsgdrD_corrvalues'].squeeze()
        SWHvalues = mat['SWHinstrcorr_J3WHALESsgdrD_swhvalues'].squeeze()
        
        if SWH_ALES < np.min(SWHvalues):
            instr_corr_SWH = 0
        elif SWH_ALES > np.max(SWHvalues):          
            instr_corr_SWH = corrvalues[-1]                
        else:
            f = interpolate.interp1d(SWHvalues, corrvalues, kind='cubic')
            instr_corr_SWH = f(SWH_ALES)  
    
    return instr_corr_SWH, f  
