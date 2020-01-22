# -*- coding: utf-8 -*-
"""
Created on November 2018

The code launches the WHALES retracker using original mission files

Works for the following missions: 

@author: Marcello Passaro
"""

#import cmath
import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use("Agg")
#import matplotlib.pyplot as plt
import scipy.io
import os
#import matplotlib.animation as manimation
import time
from compute_instr_corr_SWH_WHALES import compute_instr_corr_SWH_WHALES
#import sys
# from read_functions import wf_reader

from Retracker_MP import *

from WHALES_withRangeAndEpoch import *

from scipy.io import matlab




## EXPERIMENTS: '/DGFI34/data/ALES_Python/Reprocessing/prove/'
saving_directory= '' 

filename='JA3_GPS_2PdP054_149_20170801_175030_20170801_184643.nc'

mission='jason3'

cal2='on'

add_instr_corr_SWH = 'yes'
import_weights= 'yes'

interpolation_factor=1 #Choose whether waveform is oversampled or not


J2_filter = np.loadtxt('J2_MeanFilterKu')
J2_filter_norm = J2_filter / np.mean(J2_filter[11:115])  
J3_filter = np.loadtxt('J3_MeanFilterKu')
J3_filter_norm = J3_filter / np.mean(J3_filter[11:115])       
         
if add_instr_corr_SWH == 'yes' :
    my_path_instr_corr_SWH='SWHinstrcorr_WHALES_jason3SGDRd.mat' 
    
if import_weights == 'yes' :
    my_path_weights='weights.mat'     
    mat_weights = matlab.loadmat(my_path_weights)
    residual_std=np.squeeze(mat_weights['residual_tot'])
    flag_edges=np.squeeze(mat_weights['flag_edges'])
     


# 2) FUNCTION DEFINITIONS

def moving_average(a, n) :
    a = np.concatenate( ( np.zeros( (n+1,) )  , a , np.zeros( (n,) )    )   );
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def wf_reader(filename):

 
        
       
    S= netCDF4.Dataset(filename,'r')
        

           

    return S  


####### ---------------------------------------------------------------######



# 3) Launcher
counting=0 # Only needed to check which is the first waveform to be reprocessed, in order to use the SSB interpolation only once
counting_swh=0 # Only needed to check which is the first waveform to be reprocessed, in order to use the swh_instr_corr interpolation only once



S  = wf_reader(filename)

#HERE ADAPT FOR JASON-3 NAMING
S_time=np.ma.getdata( S.variables['time_20hz'][:] )
S_height=np.ma.getdata( S.variables['alt_20hz'][:] )
S_swh=np.ma.getdata( S.variables['swh_20hz_ku'][:] )
S_tracker=np.ma.getdata( S.variables['tracker_20hz_ku'][:] )
S_range=np.ma.getdata( S.variables['range_20hz_ku'][:] )
S_waveform=np.ma.getdata( S.variables['waveforms_20hz_ku'][:] )
S_lat=np.ma.getdata( S.variables['lat_20hz'][:] )
S_lon=np.ma.getdata( S.variables['lon_20hz'][:] )
S_landmask=np.ma.getdata( S.variables['surface_type'][:] )
S_offnadir=np.ma.getdata( S.variables['off_nadir_angle_wf_20hz_ku'][:] )


S_atmos_corr=np.ma.getdata( S.variables['atmos_corr_sig0_ku'][:] )
#This field is at 1-Hz, so it has to be reshaped
S_atmos_corr=np.transpose(np.tile(S_atmos_corr,(np.shape(S_time)[1],1)))

S_scaling_factor=np.ma.getdata( S.variables['scaling_factor_20hz_ku'][:] )
 


# WHALES RETRACKING ATTEMPT 
landmask=np.empty(np.shape(S_time))*np.nan

swh_WHALES=np.empty(np.shape(S_time))*np.nan

Err_WHALES=np.empty(np.shape(S_time))*np.nan
Epoch_WHALES=np.empty(np.shape(S_time))*np.nan
Amplitude_WHALES=np.empty(np.shape(S_time))*np.nan

sigma0_WHALES=np.empty(np.shape(S_time))*np.nan


time_20hz = np.empty(np.shape(S_time))*np.nan

altitude   =np.empty(np.shape(S_time))*np.nan           
range_WHALES=np.empty(np.shape(S_time))*np.nan


swh_WHALES_instr_corr=np.empty(np.shape(S_time))*np.nan


for index_waveforms_row in np.arange(0,np.shape(S_time)[0],1):  #np.arange(0,np.shape(S_time)[0],1)
    
    print index_waveforms_row    
    
    for index_waveforms_col in np.arange(0,np.shape(S_time)[1],1):

        print index_waveforms_col

        #landmask[index_waveforms] = parameters_landmask['mask'][index_waveforms]            
        
        print "Retracking waveform  "+ str(index_waveforms_row) +  "  of  " + str(np.size(S_time))  
        #str(cycle_index) +  "  of  " + str(np.size(cycle_vector)) +"...pass...  "+ str(path_index) +  "  of  " + str(np.size(path_vector))
        
        
         
        if mission=='jason3' or mission=='jason2':
            
            antenna_ref_point_correction=0.18092 # See for example: http://www.aviso.oceanobs.com/fileadmin/documents/data/tools/JA2_GDR_D_release_note.pdf
                                          # The full reference is  Desjonquerès, J. D., and N. Picot. "OSTM/JASON-2 absolute bias technical note." CNES internal document TP3-JPOS3-NT-1627-CNES (2011).
    
                   
            input = {}
            if cal2=='on'  :                      
                if mission=='jason3' :
                    ' waveform '
                    input['waveform'] = S_waveform[index_waveforms_row,index_waveforms_col,:]/J3_filter_norm[11:115]       
                elif mission=='jason2':
                    ' waveform '
                    input['waveform'] = S_waveform[index_waveforms_row,index_waveforms_col,:]/J2_filter_norm[11:115] 
            else:
                input['waveform'] = S_waveform[index_waveforms_row,index_waveforms_col,:]
                
                    
    
            ' reference bin '
            input['refbin'] =  52.0
            ' bin size in [ns] '
            input['binsize'] = 2e-09
            ' raw range in [m] '
            input['uralt'] = S_tracker[index_waveforms_row,index_waveforms_col]
            ' Doppler correction for range '
            input['doppler'] = 0 #Not applied here
            ' hsat '        
            input['hsat'] = S_height[index_waveforms_row,index_waveforms_col]
            ' mission '
            if mission=='jason3' :
                input['mission'] = 'jason3'
            elif mission=='jason2':
                input['mission'] = 'jason2'
            ' beamwidth in degree '
            
            input['theta0'] = 1.29 
            ' off nadir angle in degree ' 
            input['xi'] = S_offnadir[index_waveforms_row,index_waveforms_col]
            
            if import_weights == 'yes' :     
                input['weights_flag']=flag_edges
                input['weights']=residual_std
                                        
                
            retracker = WHALES_withRangeAndEpoch(input)
                    
                
     
                
                

            #Quality flag of WHALES, based on the normalised fitting error on the leading edge
            if (retracker.Error) > 0.3 and (np.isnan(retracker.Error)==0):                        
                Err_WHALES[index_waveforms_row,index_waveforms_col] = 1
            elif retracker.Error<=0.3 :
                Err_WHALES[index_waveforms_row,index_waveforms_col] = 0
    
            swh_WHALES[index_waveforms_row,index_waveforms_col]               =retracker.SWH
            
            Epoch_WHALES[index_waveforms_row,index_waveforms_col]               =retracker.Epoch
    
            if mission in ['envisat','envisat_over']:
                sigma0_WHALES[index_waveforms_row,index_waveforms_col]            =retracker.Amplitude+S_atmos_corr[index_waveforms_row,index_waveforms_col]+ S_scaling_factor[index_waveforms_row,index_waveforms_col] -33.1133
            elif mission in ['jason2','jason1','saral','saral_igdr','jason3']:
                sigma0_WHALES[index_waveforms_row,index_waveforms_col]            =retracker.Amplitude+S_atmos_corr[index_waveforms_row,index_waveforms_col]+ S_scaling_factor[index_waveforms_row,index_waveforms_col]
            
            range_WHALES[index_waveforms_row,index_waveforms_col]             =retracker.range
            
            Amplitude_WHALES[index_waveforms_row,index_waveforms_col] = retracker.Norm_Amplitude
    
        
                
     
    
            # APPLICATION OF INSTRUMENTAL CORRECTION FOR SWH           
            if add_instr_corr_SWH == 'yes' :
                if counting_swh == 0: 
                    swh_WHALES_instr_corr[index_waveforms_row,index_waveforms_col] ,interpolator_instr_corr_SWH=compute_instr_corr_SWH_WHALES(swh_WHALES[index_waveforms_row,index_waveforms_col],my_path_instr_corr_SWH,mission)
                    counting_swh = 1
                else :
                    swh_WHALES_instr_corr[index_waveforms_row,index_waveforms_col] ,interpolator_instr_corr_SWH=compute_instr_corr_SWH_WHALES(swh_WHALES[index_waveforms_row,index_waveforms_col],my_path_instr_corr_SWH,mission,interpolator_instr_corr_SWH)               




#NETCDF CODE
w_nc_fid = Dataset(saving_directory+'try.nc', 'w', format='NETCDF3_CLASSIC')              
w_nc_fid.createDimension('time', np.shape(time_20hz)[0])
w_nc_fid.createDimension('records', np.shape(time_20hz)[1])

w_nc_var = w_nc_fid.createVariable('time_20hz', 'f8', ('time','records'),zlib=True)
w_nc_var.setncatts({'long_name': u"time_20hz",\
                    'units': u"s",\
                    'comment': u"time in seconds"})
w_nc_fid.variables['time_20hz'][:] = S_time           

w_nc_var = w_nc_fid.createVariable('lat_20hz', 'f8', ('time','records'),zlib=True)
w_nc_var.setncatts({'long_name': u"lat_20hz",\
                    'units': u"deg",\
                    'comment': u" "})
w_nc_fid.variables['lat_20hz'][:] = S_lat   

w_nc_var = w_nc_fid.createVariable('lon_20hz', 'f8', ('time','records'),zlib=True)
w_nc_var.setncatts({'long_name': u"lon_20hz",\
                    'units': u"deg",\
                    'comment': u" "})
w_nc_fid.variables['lon_20hz'][:] = S_lon   

w_nc_var = w_nc_fid.createVariable('swh_WHALES_20hz', 'f8', ('time','records'),zlib=True)
w_nc_var.setncatts({'long_name': u"swh_WHALES_20hz",\
                    'units': u"m",\
                    'comment': u" "})
w_nc_fid.variables['swh_WHALES_20hz'][:] = swh_WHALES

w_nc_var = w_nc_fid.createVariable('swh_WHALES_instr_corr_20hz', 'f8', ('time','records'),zlib=True)
w_nc_var.setncatts({'long_name': u"swh_WHALES_instr_corr_20hz",\
                    'units': u"m",\
                    'comment': u" "})
w_nc_fid.variables['swh_WHALES_instr_corr_20hz'][:] = swh_WHALES_instr_corr

w_nc_var = w_nc_fid.createVariable('sigma0_WHALES_20hz', 'f8', ('time','records'),zlib=True)
w_nc_var.setncatts({'long_name': u"sigma0_WHALES_20hz",\
                    'units': u"dB",\
                    'comment': u" "})
w_nc_fid.variables['sigma0_WHALES_20hz'][:] = sigma0_WHALES
  
w_nc_var = w_nc_fid.createVariable('range_WHALES_20hz', 'f8', ('time','records'),zlib=True)
w_nc_var.setncatts({'long_name': u"range_WHALES_20hz",\
                    'units': u"m",\
                    'comment': u" "})
w_nc_fid.variables['range_WHALES_20hz'][:] = range_WHALES

w_nc_var = w_nc_fid.createVariable('epoch_WHALES_20hz', 'f8', ('time','records'),zlib=True)
w_nc_var.setncatts({'long_name': u"epoch_WHALES_20hz",\
                    'units': u"m",\
                    'comment': u" "})
w_nc_fid.variables['epoch_WHALES_20hz'][:] = Epoch_WHALES


#
#w_nc_var = w_nc_fid.createVariable('lat_1hz', 'f8', ('time'),zlib=True,least_significant_digit=4)
#w_nc_var.setncatts({'long_name': u"Latitude",\
#                    'units': u"deg North"})
#w_nc_fid.variables['lat_1hz'][:] = lat

w_nc_fid.close()  # close the new file




