# -*- coding: utf-8 -*-
"""
Created on Sat Dec. 21 2024

@author: M. Passaro & F. ardhuin
"""

"""altimeters_parameters.py: A Python module for LRM waveform models
"""
#================================================================
# Imports
#----------------------------------------------------------------
import numpy as np

######################  Defines parameters
def  instrument_parameters(mission)  :
# Default values for the Jason series 
    tau=3.125 #gate spacing in ns
    SigmaP=0.513*tau
    nump=90
    total_gate_number=104
    nominal_tracking_gate=31     
                
# Mission-dependent parameters and files to be loaded
    if mission.lower() == 'ers2_r' or mission.lower() == 'ers2_r_2cm':           
        tau=3.03
        Theta=1.3 *np.pi/180
        SigmaP=0.513*tau    
        nump=50
    if mission in ['envisat']:
        Theta=1.35 *np.pi/180 #modified http://www.aviso.oceanobs.com/fileadmin/documents/OSTST/2010/oral/PThibaut_Jason2.pdf  % The antenna 3dB bandwidth (degrees transformed in radians)
        SigmaP=0.53*tau #from Gomez Enri 2006. Otherwise use:%1.6562; %ns =0.53*3.125ns
        total_gate_number=128
        nominal_tracking_gate=45
    elif mission in ['jason1']:
        Theta=1.29 *np.pi/180
    elif mission in ['jason2']:
        Theta=1.29 *np.pi/180
    elif mission in ['jason3', 'jason3f','jason3f2','swot']:
        Theta=1.29 *np.pi/180
    elif mission.lower() in ['altika', 'saral', 'saral_igdr']:
        tau=3.125*320/480      
        Theta=0.605 *np.pi/180
        SigmaP=0.513*tau
        nump=96  # Steunou et al. 2015 
        total_gate_number=128
        nominal_tracking_gate=51          
    elif mission in ['cs2_lrm']:
        tau=3.125 #gate spacing in ns
        Theta=1.1992 *np.pi/180 #modified http://www.aviso.oceanobs.com/fileadmin/documents/OSTST/2010/oral/PThibaut_Jason2.pdf  % The antenna 3dB bandwidth (degrees transformed in radians)
        SigmaP=0.513*tau    
        nump=95 
        total_gate_number=128                
        nominal_tracking_gate=64          
    return Theta,tau,SigmaP,nump,total_gate_number,nominal_tracking_gate

def  setpaths_corrections(mission)  :
    if mission in ['envisat']:
        my_path_instr_corr_SWH = ''
        my_path_weights = 'weights/weights_n1.mat'
    elif mission in ['jason1']:
        my_path_instr_corr_SWH = 'instr_corr/SWHinstrcorr_MLE4_jason1SGDRc.mat'
        my_path_weights = 'weights/weights.mat'
    elif mission in ['jason2']:
        my_path_instr_corr_SWH = 'instr_corr/SWHinstrcorr_WHALES_jason2SGDRd.mat'
        my_path_weights = 'weights/weights_J2.pkl'
    elif mission in ['jason3', 'jason3f','jason3f2','swot']:
        my_path_instr_corr_SWH = 'instr_corr/SWHinstrcorr_WHALES_jason3SGDRd.mat'
        my_path_weights = 'weights/weights_J2.pkl'
    elif mission.lower() in ['altika', 'saral', 'saral_igdr']:
        my_path_instr_corr_SWH = ''
        my_path_weights = 'weights/weights_alt.mat'
    elif mission in ['cs2_lrm']:
        my_path_instr_corr_SWH = ''
        my_path_weights = 'weights/weights_cs2_lrm.mat'
    return my_path_instr_corr_SWH,my_path_weights

################################################################################################################    
def  processing_choices(mission)  :
    thra=0.06   # threshold for normalized waveform at second range gate of leading edge
    # This second threshold was introduced by FA (to recover previous behavior, set thrb to 0) . 
    thrb=0.8   # threshold for lowest normalized waveform value beyond which the leading edge may stop. 
    minHs=0.2
    maxHs=30
    if mission.lower() == 'envisat':
                noisegates=np.arange(4,10); #gates used to estimate Thermal Noise
                startgate=4                 #First gate to be considered in the retracking window
                ALEScoeff0=2.45             #experimental values for SWH. it is the constant term in the definition of the number of gates to be considered in the retracking
                                #after the middle of the leading edge
                ALEScoeff1=4.05 #This is the slope of the WHALES relationship between tolerance of precision and width of the subwaveform   
                Err_tolerance_vector=0.3; #Tolerance on the (normalised) fitting error of the waveform. It can be used, for example,
                                                        #to retrack the same waveform in a different way if fitting performances are not satisfactory

    elif mission.lower() == 'saral' or mission.lower() == 'saral_igdr':
                noisegates=10+np.arange(4,10); #gates used to estimate Thermal Noise # changed by Marine De Carlo
                startgate=4 #First gate to be considered in the retracking window
                ALEScoeff0=2.94 #experimental values for SWH. it is the constant term in the definition of the number of gates to be considered in the retracking
                                #after the middle of the leading edge
                ALEScoeff1=3.56 #This is the slope of the WHALES relationship between tolerance of precision and width of the subwaveform   
                Err_tolerance_vector=0.3; #Tolerance on the (normalised) fitting error of the waveform. It can be used, for example,
                                                        #to retrack the same waveform in a different way if fitting performances are not satisfactory              

   
    elif mission.lower() == 'jason2' or mission.lower() == 'jason1' or mission.lower() == 'jason3' or mission.lower() == 'swot': 
                noisegates=np.arange(0,6); #gates used to estimate Thermal Noise
                startgate=1 #First gate to be considered in the retracking window
                ALEScoeff0=3.89 #experimental values for SWH. it is the constant term in the definition of the number of gates to be considered in the retracking
                                #after the middle of the leading edge
                ALEScoeff1=3.86 #This is the slope of the WHALES relationship between tolerance of precision and width of the subwaveform   
                Err_tolerance_vector=0.3; #Tolerance on the (normalised) fitting error of the waveform. It can be used, for example,
                                                        #to retrack the same waveform in a different way if fitting performances are not satisfactory

    elif mission.lower() == 'cs2_lrm' :
                noisegates=np.arange(4,10); #gates used to estimate Thermal Noise
                startgate=4 #First gate to be considered in the retracking window
                ALEScoeff0=3.68 #experimental values for SWH. it is the constant term in the definition of the number of gates to be considered in the retracking
                                #after the middle of the leading edge
                ALEScoeff1=3.36 #This is the slope of the WHALES relationship between tolerance of precision and width of the subwaveform   
                Err_tolerance_vector=0.3; #Tolerance on the (normalised) fitting error of the waveform. It can be used, for example,
                                                        #to retrack the same waveform in a different way if fitting performances are not satisfactory                 


    elif mission.lower() == 'ers2_r': 
                print("Mission not yet supported")
                sys.exit(0)
                
    elif mission.lower() == 'ers2_r_2cm':
                print("Mission not yet supported")
                sys.exit(0)

    else:
                print("unknown mission")
                sys.exit(0)

    return(noisegates,startgate,ALEScoeff0,ALEScoeff1,Err_tolerance_vector,thra,thrb,minHs,maxHs) 

