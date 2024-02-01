#!/bin/python

import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import traceback

from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.integrate import quad

#from scipy.integrate import simps

#from Database import *
#from Functions import *

class Retracker_MP:
    
    debug = 1
        
    ' Retracker configuration '
    waveform = None
    refbin = None
    binsize = None
    uralt = None
    threshold = None
    scale = 0
    mission = None
    hsat = None
    theta0 = None
    ref_range = None
    option = None
    factor = None
    #weights = None #Weights used for the estimation of SWH
    #weights_flag = None #Flag that identifies start and end of the leading edge in the weight distribution
    
    ' Retracker results '
    model_x_m = None
    model_x = None
    model_y = None
    
    model = None
    range = None
    leading_edge = None    
    error = None
    b1 = None
    b2 = None
    b3 = None
    b4 = None
    b5 = None
    
    # ADDITIONAL RETRACKER OUTPUT FROM ALES    
    gate1 = None #MP added initial gate of the leading edge
    gate2 = None #MP added end gate of the leading edge
    Wt_all_yang = None  #MP added normalised fitted waveform - first pass
    Wt_all_LESfive = None #MP added normalised fitted waveform - second pass  
    D = None #MP added Normalised real waveform
    Epoch = None #MP added Epoch (m)
    SWH = None #MP added SWH (m)
    Amplitude = None  #MP added normalised amplitude
    Norm_Amplitude = None
    uncssh = None     #MP added Orbit Altitude - Retracked Range
    classifier = None
    interpolator = None      
    Sigma= None #MP added Rising Time of the leading edge
    
    ' speed of light [m/s] '
    c = 299792458
    
    def __init__(self,config):        
        if config.has_key('waveform'):
            self.waveform = config['waveform']
            for i in range(0,len(self.waveform)):
                self.waveform[i] = float(self.waveform[i])
                
        if config.has_key('waveform_last'):
            self.waveform_last = config['waveform_last']
            for i in range(0,len(self.waveform_last)):
                self.waveform_last[i] = float(self.waveform_last[i])
        if config.has_key('waveform_current'):
            self.waveform_current = config['waveform_current']
            for i in range(0,len(self.waveform_current)):
                self.waveform_current[i] = float(self.waveform_current[i])                
        if config.has_key('waveform_next'):
            self.waveform_next = config['waveform_next']
            for i in range(0,len(self.waveform_next)):
                self.waveform_next[i] = float(self.waveform_next[i])
                
        if config.has_key('refbin'):
            self.refbin = config['refbin']
        if config.has_key('binsize'):
            self.binsize = config['binsize']
        if config.has_key('uralt'):
            self.uralt = config['uralt']
        if config.has_key('threshold'):
            if float(config['threshold']) < 0.0 or float(config['threshold']) > 100.0:
                print("Error: Threshold must be between 0 and 100%")
                sys.exit(0)
            self.threshold = float(config['threshold']) / 100.0
        if config.has_key('scale') and config['scale'] == 1:
            self.scale = 1
        
        ' improved threshold '
        if config.has_key('ref_range'):
            self.ref_range = config['ref_range']
        if config.has_key('option'):
            self.option = config['option']
        
        if config.has_key('factor'):
            self.factor = config['factor']
                
        #MP added
        ' ales_rev'
        if config.has_key('mission'):
            self.mission = config['mission']
        if config.has_key('hsat'):
            self.hsat = config['hsat']
        if config.has_key('theta0'):
            self.theta0 = config['theta0']
        if config.has_key('xi'):
            self.xi = config['xi'] 
        if config.has_key('doppler'):
            self.doppler = config['doppler'] 
        if config.has_key('Epoch'):
            self.Epoch = config['Epoch']
        if config.has_key('SWH'):
            self.SWH = config['SWH']
        if config.has_key('Amplitude_initial'):
            self.Amplitude_initial = config['Amplitude_initial']        
        if config.has_key('classifier'):
            self.classifier = config['classifier']  
        if config.has_key('weights'):
            self.weights = config['weights']
        if config.has_key('weights_flag'):
            self.weights_flag = config['weights_flag']            
            
        #end MP addition    
        
        ' debug flag '
        if config.has_key('debug'):
            self.debug = int(config['debug'])            
            
        if self.__class__.__name__ == "Threshold" and self.threshold == None:
            print("Error: Threshold retracker needs argument 'threshold' [0, ..., 100]")
            sys.exit(0)
        if self.__class__.__name__ == "ImprovedThreshold" and self.threshold == None:
            print("Error: ImprovedThreshold retracker needs argument 'threshold' [0, ..., 100]")
            sys.exit(0)
            
    def plot(self):
                
        x_data = []
        y_data = []        
        for i in range(0,len(self.waveform)):
            x_data.append(i)
            y_data.append(self.waveform[i])
        plt.plot(x_data, y_data, 'b.-')
        
        x_model = []
        y_model = []
        if self.model != None:
            for i in range(0,len(self.model)):
                x_model.append(i)
                y_model.append(self.model[i])
        plt.plot(x_model, y_model, 'r.-')    
                
        if self.model_x != None and self.model_y != None:            
            plt.plot(self.model_x, self.model_y, 'r-')    
    
    
        x_leading_edge = []
        y_leading_edge = []
        if self.leading_edge != None:
            x_leading_edge.append(self.leading_edge)
            y_leading_edge.append(0)

            x_leading_edge.append(self.leading_edge)
            y_leading_edge.append(max(self.waveform))
        plt.plot(x_leading_edge, y_leading_edge, 'g.-')    
        
        plt.ylabel('Power')
        plt.xlabel('Bin')

        plt.show()
    
    def scale_waveform(self, waveform):
            
        max_val = max(waveform)
        if max_val == 0:
            return waveform
        for i in range(0,len(waveform)):
            waveform[i] = waveform[i]/max_val
        return waveform
        
    def stat_waveform(self):        
        info = {}
        info['max_power'] = None
        info['pos_max_power'] = None
        
        if self.waveform != None:
            info['max_power'] = max(self.waveform)
            for i in range (0,len(self.waveform)):
                if self.waveform[i] == info['max_power']:
                    info['pos_max_power'] = i
                    break
        return info
    
    def get_model(self):
        print("No model available for this retracker")
    
