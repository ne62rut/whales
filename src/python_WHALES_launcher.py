# -*- coding: utf-8 -*-
"""
Created on November 2018

The code launches the WHALES retracker using original mission files

Works for the following missions: Jason, Jason-2, Jason-3, SARAL, Cryosat 2 (LRM),  ... 

@author: Marcello Passaro

Examples of use from command line: 
python python_WHALES_launcher.py -m jason3f -i JA3_GPS_2PfP342_001_20230609_173418_20230609_183031.nc -o WHALES
python python_WHALES_launcher.py -m swot   -i SWOT_GPS_2PfP549_004_20230611_125241_20230611_134346.nc -w 2 -o WHALES2
python python_WHALES_launcher.py -m saral -s 3 -i SRL_GPS_2PfP001_0641_20130405_141055_20130405_150113.CNES.nc -o WHALES
python python_WHALES_launcher.py -m swot -s 3 -w 2  -i SWOT_GPS_2PfP012_567_20240326_230842_20240327_000009.nc -o WHALES

        2024-07-19: added possibility to use 1/waveform for the weights               : F. Ardhuin
	2024-12-20: adds optional smoothing (gives better results for SARAL and SWOT) : M. De Carlo, F. Ardhuin
"""

import argparse
import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib

matplotlib.use("Agg")
import scipy.io
import os
import time
from compute_instr_corr_SWH_WHALES import compute_instr_corr_SWH_WHALES

from Retracker_MP import *

from WHALES_withRangeAndEpoch import *

from import_weights_npz import *

from altimeters_parameters import *


def get_options():
    parser = argparse.ArgumentParser(
        description='Retrack a SGDR file with WHALES')

    parser.add_argument(
        '-m', '--mission', type=str,
        choices=['ers1', 'ers2', 'envisat', 'jason1', 'jason2', 'jason3', 'saral', 'cs2_lrm',
                 'jason3f', 'jason3f2','cfosat','swot','sentinel6_lrm'],
        help='satellite mission'
    )
    parser.add_argument(
        '-i', '--input', type=str,
        help='path to the SGDR file'
    )
    parser.add_argument(
        '-o', '--output', type=str, default='.',
        help='path to the output repository'
    )
    parser.add_argument(
        '-w', '--weights', type=int, default=1,
        choices=[0,1,2,3],
        help='weights to be used for retracking: 0=constant, 1=input file, 2=1/brown, 3=1/brown**2'
    )
    parser.add_argument(
        '-c', '--costfunction', type=str, default='LS',
        help='cost function for retracking: LS=least squares, ML=maximum likelihood'
    )
    parser.add_argument(
        '-W', '--weightoutofsub', type=float, default=1.0,
        help='value of 1/weight outside of subwaveform: default = 1.0'
    )
    parser.add_argument(
        '-d', '--debug', type=str, default='0',
        help='-d=1 will output waveforms and weights'
    )
    parser.add_argument(
        '-s', '--smooth', type=int, default=0,
        help='size of kernel for smoothing waveform before leading edge detection: for SARAL use 3, for ERS use 5.'
    )
    parser.add_argument(
        '-S', '--Smooth', type=int, default=0,
        help='size of kernel for smoothing waveform before leading edge detection: this is for SWH-dependent smoothing applied for SWH > 8m'
    )

    return parser.parse_args()


options = get_options()
filename = options.input
mission = options.mission
weights_type = options.weights
costfunction = options.costfunction
weight_outsub = options.weightoutofsub
debug=options.debug
smooth=options.smooth
smooth_above_8m=options.Smooth
smooth_SWH_val=8
smooth_SWH_numpoints=20

print('weight type:',weights_type)
print('constfunction:',costfunction)
print('weight out of sub:',weight_outsub)
saving_directory = options.output

saving_name = os.path.join(saving_directory, os.path.basename(filename))

cal2 = 'on'

add_instr_corr_SWH = 'no'
import_weights = 'yes'

Theta,tau,SigmaP,nump,total_gate_number,nominal_tracking_gate=instrument_parameters(mission)
my_path_instr_corr_SWH,my_path_weights=setpaths_corrections(mission) 
noisegates,startgate,ALEScoeff0,ALEScoeff1,Err_tolerance_vector,thra,thrb,minHs,maxHs,noisemin,interpolation_factor,index_originalbins,gate_number_processing=processing_choices(mission)

if my_path_instr_corr_SWH != '':
    my_path_instr_corr_SWH = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), my_path_instr_corr_SWH)
if my_path_weights != '':
    my_path_weights = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), my_path_weights)

if import_weights == 'yes':
    l1=len(my_path_weights)
    ext=my_path_weights[l1-3:l1]
    if (ext=='pkl'):
         SWH_values,residual_std,flag_edges=import_weights_pkl(my_path_weights)
    elif (ext=='npz'):
         SWH_values,residual_std,flag_edges=import_weights_npz(my_path_weights)
    else:
         residual_std,flag_edges=import_weights_mat(my_path_weights)


# 2) FUNCTION DEFINITIONS

def moving_average(a, n):
    a = np.concatenate((np.zeros((n + 1,)), a, np.zeros((n,))));
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def wf_reader(filename):
    S = netCDF4.Dataset(filename, 'r')

    return S


####### ---------------------------------------------------------------######


# 3) Launcher
counting_swh = 0  # Only needed to check which is the first waveform to be reprocessed, in order to use the swh_instr_corr interpolation only once

S = wf_reader(filename)

if mission in ['jason1', 'jason2', 'jason3']:
    # Getting 20 Hz data
    S_time = np.ma.getdata(S.variables['time_20hz'][:])
    S_height = np.ma.getdata(S.variables['alt_20hz'][:])
    S_swh = np.ma.getdata(S.variables['swh_20hz_ku'][:])
    S_tracker = np.ma.getdata(S.variables['tracker_20hz_ku'][:])
    S_range = np.ma.getdata(S.variables['range_20hz_ku'][:])
    S_waveform = np.ma.getdata(S.variables['waveforms_20hz_ku'][:])
    S_lat = np.ma.getdata(S.variables['lat_20hz'][:])
    S_lon = np.ma.getdata(S.variables['lon_20hz'][:])
    S_landmask = np.ma.getdata(S.variables['surface_type'][:])
    if 'off_nadir_angle_wf_20hz_ku' in S.variables:
        S_offnadir = np.ma.getdata(S.variables['off_nadir_angle_wf_20hz_ku'][:])
    else:
        S_offnadir = np.ma.getdata(
            S.variables['off_nadir_angle_wf_ku'][:])
        S_offnadir = np.transpose(np.tile(S_offnadir, (np.shape(S_time)[1], 1)))

    S_atmos_corr = np.ma.getdata(S.variables['atmos_corr_sig0_ku'][:])
    # This field is at 1-Hz, so it has to be reshaped
    S_atmos_corr = np.transpose(np.tile(S_atmos_corr, (np.shape(S_time)[1], 1)))

    S_scaling_factor = np.ma.getdata(S.variables['scaling_factor_20hz_ku'][:])

elif mission in ['jason3f']:
    # AVISO SGDR version F
    S_time = np.ma.getdata(S['data_20'].variables['time'][:])
    S_time = np.reshape(S_time, (np.shape(S_time)[0], 1))

    S_height = np.ma.getdata(S['data_20'].variables['altitude'][:])
    S_height = np.reshape(S_height, (np.shape(S_time)[0], 1))

    S_swh = np.ma.getdata(S['data_20']['ku'].variables['swh_ocean'][:])
    S_swh = np.reshape(S_swh, (np.shape(S_time)[0], 1))

    S_tracker = np.ma.getdata(
        S['data_20']['ku'].variables['tracker_range_calibrated'][:])
    S_tracker = np.reshape(S_tracker, (np.shape(S_time)[0], 1))

    S_range = np.ma.getdata(S['data_20']['ku'].variables['range_ocean'][:])
    S_range = np.reshape(S_range, (np.shape(S_time)[0], 1))

    S_waveform = np.ma.getdata(
        S['data_20']['ku'].variables['power_waveform'][:])
    S_waveform = np.reshape(
        S_waveform, (np.shape(S_time)[0], 1, np.shape(S_waveform)[1]))
    print(S_waveform.shape)

    S_lat = np.ma.getdata(S['data_20'].variables['latitude'][:])
    S_lat = np.reshape(S_lat, (np.shape(S_time)[0], 1))

    S_lon = np.ma.getdata(S['data_20'].variables['longitude'][:])
    S_lon = np.reshape(S_lon, (np.shape(S_time)[0], 1))

    S_landmask = np.ma.getdata(
        S['data_20'].variables['surface_classification_flag'][:])
    S_landmask = np.reshape(S_landmask, (np.shape(S_time)[0], 1))

    S_offnadir = np.ma.getdata(
        S['data_20']['ku'].variables['off_nadir_angle_wf_ocean'][:])
    S_offnadir = np.reshape(S_offnadir, (np.shape(S_time)[0], 1))

    atmos_corr = np.ma.getdata(
        S['data_01']['ku'].variables['sig0_cor_atm'][:])
    # This field is at 1-Hz, so it has to be reshaped
    S_atmos_corr = atmos_corr[S['data_20'].variables['index_1hz_measurement']]
    S_atmos_corr = np.reshape(S_atmos_corr, (np.shape(S_time)[0], 1))

    S_scaling_factor = np.ma.getdata(
        S['data_20']['ku'].variables['sig0_scaling_factor'][:])
    S_scaling_factor = np.reshape(S_scaling_factor, (np.shape(S_time)[0], 1))

    mission = 'jason3'

elif mission in ['jason3f2','swot']:
    # AVISO SGDR version F: need to check with actual j3 file 
    S_tim1 = np.ma.getdata(S['data_01'].variables['time'][:])
    ind20  = np.ma.getdata(S['data_01'].variables['index_first_20hz_measurement'][:])
    num20  = np.ma.getdata(S['data_01'].variables['numtotal_20hz_measurement'][:])
    S_time1   = np.ma.getdata(S['data_20'].variables['time'][:])
    S_height1 = np.ma.getdata(S['data_20'].variables['altitude'][:])
    S_swh1    = np.ma.getdata(S['data_20']['ku'].variables['swh_ocean'][:])
    S_tracker1= np.ma.getdata(S['data_20']['ku'].variables['tracker_range_calibrated'][:])
    S_range1= np.ma.getdata(S['data_20']['ku'].variables['range_ocean'][:])
    S_waveform1= np.ma.getdata(S['data_20']['ku'].variables['power_waveform'][:])
    S_lat1 = np.ma.getdata(S['data_20'].variables['latitude'][:])
    S_lon1 = np.ma.getdata(S['data_20'].variables['longitude'][:])
    S_landmask1 = np.ma.getdata(
        S['data_20'].variables['surface_classification_flag'][:])
    S_offnadir1 = np.ma.getdata(
        S['data_20']['ku'].variables['off_nadir_angle_wf_ocean'][:])
    S_atmos_corr1 = np.ma.getdata(
        S['data_01']['ku'].variables['sig0_cor_atm'][:])
    S_scaling_factor1 = np.ma.getdata(
        S['data_20']['ku'].variables['sig0_scaling_factor'][:])
        
    n01=np.shape(S_tim1)[0]
    [n20,nr]=np.shape(S_waveform1)
    print('shape:',n01)
    S_time=np.zeros((n01,20))
    S_height=np.zeros((n01,20))
    S_swh=np.zeros((n01,20))
    S_tracker=np.zeros((n01,20))
    S_range=np.zeros((n01,20))
    S_waveform=np.zeros((n01,20,nr))
    S_lat=np.zeros((n01,20))
    S_lon=np.zeros((n01,20))
    S_landmask=np.zeros((n01,20))
    S_offnadir=np.zeros((n01,20))
    S_atmos_corr=np.zeros((n01,20))
    S_scaling_factor=np.zeros((n01,20))
    for i01 in np.arange(np.shape(S_tim1)[0]):
        i2=ind20[i01]
        n2=num20[i01]
        S_time[i01,0:n2] = S_time1[i2:i2+n2]
        S_height[i01,0:n2] =S_height1[i2:i2+n2]
        S_swh[i01,0:n2] =S_swh1[i2:i2+n2]
        S_tracker[i01,0:n2] =S_tracker1[i2:i2+n2]
        S_range[i01,0:n2] =S_range1[i2:i2+n2]
        S_waveform[i01,0:n2,:] =S_waveform1[i2:i2+n2,:]
        S_lon[i01,0:n2] =S_lon1[i2:i2+n2]
        S_lat[i01,0:n2] =S_lat1[i2:i2+n2]
        S_landmask[i01,0:n2] =S_landmask1[i2:i2+n2]
        S_offnadir[i01,0:n2] =S_offnadir1[i2:i2+n2]
        S_atmos_corr[i01,0:n2] =S_atmos_corr1[i01]
        S_scaling_factor[i01,0:n2] =S_scaling_factor1[i2:i2+n2]

    print('Waveform array shape:',S_waveform.shape)
    if (mission=='jason3f2'): 
       mission = 'jason3'

elif mission in['sentinel6_lrm'] :

    S_group=S.groups['data_20']
    S_group=S_group['ku']
    #print(S_group.variables)
    S_time=np.ma.getdata(S_group.variables['time'][:])
    S_time=np.reshape(S_time,(np.shape(S_time)[0],1) )
    #print np.shape(S_time)
    S_height=np.ma.getdata(S_group.variables['altitude'][:])
    S_height=np.reshape(S_height,(np.shape(S_time)[0],1) )

    #S_swh=np.ma.getdata(S_group.variables['swh_ocean_20_ku'][:])

    S_tracker=np.ma.getdata(S_group.variables['tracker_range_calibrated'][:])
    S_tracker=np.reshape(S_tracker,(np.shape(S_time)[0],1) )

    #S_range=np.ma.getdata(S.variables['range_ku_l1bs_echo_sar_ku'][:])

    S_waveform_scale_factor=np.ma.getdata(S_group.variables['waveform_scale_factor'][:])
    S_waveform_scale_factor=np.reshape(S_waveform_scale_factor,(np.shape(S_time)[0],1) )

    S_waveform=np.ma.getdata(S_group.variables['power_waveform'][:])

    S_lat=np.ma.getdata(S_group.variables['latitude'][:])
    S_lat=np.reshape(S_lat,(np.shape(S_time)[0],1) )

    S_lon=np.ma.getdata(S_group.variables['longitude'][:])
    S_lon=np.reshape(S_lon,(np.shape(S_time)[0],1) )

    sig0_scaling_factor=np.ma.getdata(S_group.variables['sig0_scaling_factor'][:])
    sig0_scaling_factor=np.reshape(sig0_scaling_factor,(np.shape(S_time)[0],1) )

    ptr_main_lobe_width=np.ma.getdata(S_group.variables['ptr_main_lobe_width'][:])
    ptr_main_lobe_width=np.reshape(ptr_main_lobe_width,(np.shape(S_time)[0],1) )

    S_landmask=np.ones((len(S_tracker),1))*0 

    S_offnadir=np.ma.getdata(S_group.variables['off_nadir_roll_angle_pf'][:])
    S_offnadir=np.reshape(S_offnadir,(np.shape(S_time)[0],1) )
 

elif mission in ['ers2','ers1']:

    S_qual_wf_not_tracking=np.ma.getdata( S.variables['qual_wf_not_tracking_20hz'][:] )
    #S_time=np.reshape(S_time,(np.shape(S_time)[0],1) )

    S_time=np.ma.getdata( S.variables['time_20hz'][:] )
    #S_time=np.reshape(S_time,(np.shape(S_time)[0],1) )
    
    #Time at 1-Hz for interpolation of fields available only at 1-Hz
    S_time_1hz=np.ma.getdata( S.variables['time'][:] )
    #S_time_1hz=np.reshape(S_time_1hz,(np.shape(S_time_1hz)[0],1) )
    
    S_height=np.ma.getdata( S.variables['alt_20hz'][:] )
    #S_height=np.reshape(S_height,(np.shape(S_time)[0],1) )

    S_swh=np.ma.getdata( S.variables['swh_20hz'][:] )
    #S_swh=np.reshape(S_swh,(np.shape(S_time)[0],1) )

    S_tracker=np.ma.getdata( S.variables['tracker_range_20hz'][:] )
    #S_tracker=np.reshape(S_tracker,(np.shape(S_time)[0],1) )

    S_range=np.ma.getdata( S.variables['ocean_range_20hz'][:] )
    #S_range=np.reshape(S_range,(np.shape(S_time)[0],1) )

    S_waveform=np.ma.getdata( S.variables['ku_wf'][:] ).astype(np.float64)

    S_lat=np.ma.getdata( S.variables['lat_20hz'][:] )
    #S_lat=np.reshape(S_lat,(np.shape(S_time)[0],1) )

    S_lon=np.ma.getdata( S.variables['lon_20hz'][:] )
    #S_lon=np.reshape(S_lon,(np.shape(S_time)[0],1) )

    #S_landmask=np.ma.getdata( S.variables['surface_type'][:] )
    
    S_offnadir=np.ma.getdata( S.variables['off_nadir_angle_wf_20hz'][:] )
    #S_offnadir=np.reshape(S_offnadir,(np.shape(S_time)[0],1) )
    S_offnadir=np.zeros_like(S_offnadir) # THE OFF NADIR FIELD IN ERS2 IS EMPTY!!!!!

    S_atmos_corr=np.ma.getdata( S.variables['atmos_corr_sig0'][:] )
    #S_atmos_corr=np.reshape(S_atmos_corr,(np.shape(S_time_1hz)[0],1) )
    #This field is at 1-Hz, so it has to be reshaped
    S_atmos_corr=np.transpose(np.tile(S_atmos_corr,(np.shape(S_time)[1],1)))
    #S_atmos_corr=np.interp(S_time[:,0],S_time_1hz[:,0],S_atmos_corr[:,0])
    #S_atmos_corr=np.reshape(S_atmos_corr,(np.shape(S_time)[0],1) )

    S_scaling_factor=np.ma.getdata( S.variables['scaling_factor_20hz'][:] )
    #S_scaling_factor=np.reshape(S_scaling_factor,(np.shape(S_time)[0],1) )    

elif mission in ['envisat']:
    
    S_time=np.ma.getdata( S.variables['time_20'][:] ) 
    print(np.shape(S_time))
    S_time=np.reshape(S_time,(np.shape(S_time)[0],1) )
    print(np.shape(S_time))
    
    #Time at 1-Hz for interpolation of fields available only at 1-Hz
    S_time_1hz=np.ma.getdata( S.variables['time_01'][:] )
    S_time_1hz=np.reshape(S_time_1hz,(np.shape(S_time_1hz)[0],1) )

    S_height=np.ma.getdata( S.variables['alt_20'][:] )
    S_height=np.reshape(S_height,(np.shape(S_time)[0],1) )
    
    S_swh=np.ma.getdata( S.variables['swh_ocean_20_ku'][:] )
    S_swh=np.reshape(S_swh,(np.shape(S_time)[0],1) )
    print(np.shape(S_height))
    
    S_tracker=np.ma.getdata( S.variables['tracker_range_20_ku'][:] )
    S_tracker=np.reshape(S_tracker,(np.shape(S_time)[0],1) )
    
    S_range=np.ma.getdata( S.variables['range_ocean_20_ku'][:] )
    S_range=np.reshape(S_range,(np.shape(S_time)[0],1) )
    
    S_waveform=np.ma.getdata( S.variables['waveform_fft_20_ku'][:] )
    #S_waveform=np.reshape(S_waveform,(np.shape(S_time)[0],1) )
    
    S_lat=np.ma.getdata( S.variables['lat_20'][:] )
    S_lat=np.reshape(S_lat,(np.shape(S_time)[0],1) )
    
    S_lon=np.ma.getdata( S.variables['lon_20'][:] )
    S_lon=np.reshape(S_lon,(np.shape(S_time)[0],1) )
    
    #S_landmask=np.ma.getdata( S.variables['surf_type_20'][:] )
    #S_landmask=np.reshape(S_landmask,(np.shape(S_time)[0],1) )
    
    #OFF NADIR ANGLE FROM WAVEFORM
    #S_offnadir=np.ma.getdata( S.variables['off_nadir_angle_wf_ocean_20_ku'][:] )     #degrees^2
    #S_offnadir=np.reshape(S_offnadir,(np.shape(S_time)[0],1) )
    #S_offnadir=S_time*0.0
    
    #OFF NADIR ANGLE FROM PLATFORM
    S_offnadir=np.ma.getdata( S.variables['off_nadir_angle_pf_01'][:] )
    S_offnadir=np.reshape(S_offnadir,(np.shape(S_time_1hz)[0],1) )
    S_offnadir=np.interp(S_time[:,0],S_time_1hz[:,0],S_offnadir[:,0])
    S_offnadir=np.reshape(S_offnadir,(np.shape(S_time)[0],1) )    
    
    S_atmos_corr=np.ma.getdata( S.variables['atm_cor_sig0_01_ku'][:] )
    S_atmos_corr=np.reshape(S_atmos_corr,(np.shape(S_time_1hz)[0],1) )
    #This field is at 1-Hz, so it has to be reshaped
    S_atmos_corr=np.interp(S_time[:,0],S_time_1hz[:,0],S_atmos_corr[:,0])
    S_atmos_corr=np.reshape(S_atmos_corr,(np.shape(S_time)[0],1) )

    S_scaling_factor=np.ma.getdata( S.variables['scale_factor_20_ku'][:] )
    S_scaling_factor=np.reshape(S_scaling_factor,(np.shape(S_time)[0],1) )
    

elif mission in ['saral', 'altika']:
    S_time = np.ma.getdata(S.variables['time_40hz'][:])
    S_height = np.ma.getdata(S.variables['alt_40hz'][:])
    S_swh = np.ma.getdata(S.variables['swh_40hz'][:])
    S_tracker = np.ma.getdata(S.variables['tracker_40hz'][:])
    S_range = np.ma.getdata(S.variables['range_40hz'][:])
    S_waveform = np.ma.getdata(S.variables['waveforms_40hz'][:])
    S_lat = np.ma.getdata(S.variables['lat_40hz'][:])
    S_lon = np.ma.getdata(S.variables['lon_40hz'][:])
    S_landmask = np.ma.getdata(S.variables['surface_type'][:])

    # OFF NADIR ANGLE FROM WAVEFORM
    # S_offnadir=np.ma.getdata( S.variables['off_nadir_angle_wf_40hz'][:] )
    ##Unrealistic offnadir angles of value higher than 0.3 degrees, which would affect Range estimation (Dorandeau et al. 2004) are removed
    ##Note that in Altika there are some values put as 3267, where the retracking of the offnadir likely failed
    # index_offnadir=np.where(np.abs(S_offnadir)>0.3)[0]
    # S_offnadir[index_offnadir]=0.
    ##Off nadir angle field is filtered with an alongtrack filter of 30 seconds as suggested by Amarouche et al. (2004)

    # OFF NADIR ANGLE FROM PLATFORM
    S_offnadir = np.ma.getdata(S.variables['off_nadir_angle_pf'][:])
    # This field is at 1-Hz, so it has to be reshaped
    S_offnadir = np.transpose(np.tile(S_offnadir, (np.shape(S_time)[1], 1)))

    S_atmos_corr = np.ma.getdata(S.variables['atmos_corr_sig0'][:])
    # This field is at 1-Hz, so it has to be reshaped
    S_atmos_corr = np.transpose(np.tile(S_atmos_corr, (np.shape(S_time)[1], 1)))

    S_scaling_factor=np.ma.getdata( S.variables['scaling_factor_40hz'][:] )  
    
elif mission in ['cs2_lrm']:

    S_time=np.ma.getdata( S.variables['time_20_ku'][:] ) 
    print(np.shape(S_time))
    S_time=np.reshape(S_time,(np.shape(S_time)[0],1) )
    print(np.shape(S_time))
    
    #Time at 1-Hz for interpolation of fields available only at 1-Hz
    S_time_1hz=np.ma.getdata( S.variables['time_avg_01_ku'][:] )
    S_time_1hz=np.reshape(S_time_1hz,(np.shape(S_time_1hz)[0],1) )

    S_height = np.ma.getdata(S.variables['alt_20_ku'][:])
    S_height = np.reshape(S_height, (np.shape(S_time)[0], 1))

    # S_swh=np.ma.getdata( S.variables['swh_ocean_20_ku'][:] )
    # S_swh=np.reshape(S_swh,(np.shape(S_time)[0],1) )
    # print(np.shape(S_height))

    S_tracker = np.ma.getdata(
        S.variables['window_del_20_ku'][:] * (299792458.0) / 2.0)      # warning : light speed should be constant 
    S_tracker = np.reshape(S_tracker, (np.shape(S_time)[0], 1))

    # S_range=np.ma.getdata( S.variables['range_ocean_20_ku'][:] )
    # S_range=np.reshape(S_range,(np.shape(S_time)[0],1) )

    S_waveform = np.ma.getdata(S.variables['pwr_waveform_20_ku'][:]).astype(np.float64)
    # S_waveform=np.reshape(S_waveform,(np.shape(S_time)[0],1) )

    S_lat = np.ma.getdata(S.variables['lat_20_ku'][:])
    S_lat = np.reshape(S_lat, (np.shape(S_time)[0], 1))

    S_lon = np.ma.getdata(S.variables['lon_20_ku'][:])
    S_lon = np.reshape(S_lon, (np.shape(S_time)[0], 1))

    # S_landmask=np.ma.getdata( S.variables['surf_type_20'][:] )
    # S_landmask=np.reshape(S_landmask,(np.shape(S_time)[0],1) )

    # S_offnadir=np.ma.getdata( S.variables['off_nadir_angle_wf_ocean_20_ku'][:] )     #degrees^2
    # S_offnadir=np.reshape(S_offnadir,(np.shape(S_time)[0],1) )

    # OFF NADIR ANGLE FROM PLATFORM
    S_offnadir = np.ma.getdata(
        S.variables['off_nadir_pitch_angle_str_20_ku'][:])
    # S_offnadir=np.reshape(S_offnadir,(np.shape(S_time_1hz)[0],1) )
    # S_offnadir=np.interp(S_time[:,0],S_time_1hz[:,0],S_offnadir[:,0])
    S_offnadir = np.reshape(S_offnadir, (np.shape(S_time)[0], 1))

    # S_atmos_corr=np.ma.getdata( S.variables['atm_cor_sig0_01_ku'][:] )
    # S_atmos_corr=np.reshape(S_atmos_corr,(np.shape(S_time_1hz)[0],1) )
    ##This field is at 1-Hz, so it has to be reshaped
    # S_atmos_corr=np.interp(S_time[:,0],S_time_1hz[:,0],S_atmos_corr[:,0])
    # S_atmos_corr=np.reshape(S_atmos_corr,(np.shape(S_time)[0],1) )

    # S_scaling_factor=np.ma.getdata( S.variables['scale_factor_20_ku'][:] )
    # S_scaling_factor=np.reshape(S_scaling_factor,(np.shape(S_time)[0],1) )

nr=np.shape(S_waveform)[-1]
print(np.shape(S_waveform),'shape of S_time array:',np.shape(S_time),', nr=',nr)
# WHALES RETRACKING ATTEMPT
landmask = np.empty(np.shape(S_time)) * np.nan

swh_WHALES = np.empty(np.shape(S_time)) * np.nan
startgate_WHALES = np.empty(np.shape(S_time),dtype=np.uint16)
endgate_WHALES   = np.empty(np.shape(S_time),dtype=np.uint16)
finalgate_WHALES = np.empty(np.shape(S_time),dtype=np.uint16) 
noise_WHALES = np.empty(np.shape(S_time)) * np.nan
scale_WHALES = np.empty(np.shape(S_time)) * np.nan

Err_WHALES = np.empty(np.shape(S_time),dtype=np.uint8)
Epoch_WHALES = np.empty(np.shape(S_time)) * np.nan
Amplitude_WHALES = np.empty(np.shape(S_time)) * np.nan

sigma0_WHALES = np.empty(np.shape(S_time)) * np.nan

time_20hz = np.empty(np.shape(S_time)) * np.nan

altitude = np.empty(np.shape(S_time)) * np.nan
range_WHALES = np.empty(np.shape(S_time)) * np.nan

swh_WHALES_instr_corr = np.empty(np.shape(S_time)) * np.nan


#########################################################################################
#   Creation of NetCDF output file
#########################################################################################
w_nc_fid = Dataset(saving_name, 'w',
                   format='NETCDF3_CLASSIC')
w_nc_fid.createDimension('time', np.shape(time_20hz)[0])
w_nc_fid.createDimension('records', np.shape(time_20hz)[1])
if debug=='1':
    w_nc_fid.createDimension('gates', nr)

# global attributes: 
w_nc_fid.smooth=smooth
w_nc_fid.weight_outsub=weight_outsub
w_nc_fid.ALEScoeff0=ALEScoeff0
w_nc_fid.ALEScoeff1=ALEScoeff1
w_nc_fid.threshold_a=thra
w_nc_fid.threshold_b=thrb
w_nc_fid.maxHs=maxHs
w_nc_fid.noise_floor_is_min=noisemin

w_nc_var = w_nc_fid.createVariable('time_20hz', 'f8', ('time', 'records'),
                                   zlib=True)
w_nc_var.setncatts({'long_name': u"time_20hz", \
                    'units': u"s", \
                    'comment': u"time in seconds"})
w_nc_fid.variables['time_20hz'][:] = S_time

w_nc_var = w_nc_fid.createVariable('lat_20hz', 'f4', ('time', 'records'),
                                   zlib=True)
w_nc_var.setncatts({'long_name': u"lat_20hz", \
                    'units': u"deg", \
                    'comment': u" "})
w_nc_fid.variables['lat_20hz'][:] = S_lat

w_nc_var = w_nc_fid.createVariable('lon_20hz', 'f4', ('time', 'records'),
                                   zlib=True)
w_nc_var.setncatts({'long_name': u"lon_20hz", \
                    'units': u"deg", \
                    'comment': u" "})
w_nc_fid.variables['lon_20hz'][:] = S_lon

if debug=='1':
    w_nc_var1 = w_nc_fid.createVariable('normalized_waveform', 'f4', ('time', 'records','gates'),
                                   zlib=True)

    w_nc_var2 = w_nc_fid.createVariable('weights', 'f4', ('time', 'records','gates'),
                                   zlib=True)

    w_nc_var3 = w_nc_fid.createVariable('fitted_waveform', 'f4', ('time', 'records','gates'),
                                   zlib=True)

#
# Now looping over waveforms for retracking
# First loop is on 1 Hz data, second loop is on higher rate data (20 Hz for Jason)
#
print('size:',np.shape(S_time)[0],np.shape(S_time)[1])
for index_waveforms_row in np.arange(0,np.shape(S_time)[0], 1):
#for index_waveforms_row in np.arange(0, 10, 1):
    print("Retracking waveform group " + str(index_waveforms_row) + "  of  " +
              str(np.shape(S_time)[0]))
    for index_waveforms_col in np.arange(0, np.shape(S_time)[1], 1):
        if mission=='ers2' or mission=='ers1':
            if S_qual_wf_not_tracking[index_waveforms_row,index_waveforms_col]==1:
                print('Periodic peaky waveform, excluded using qual_wf_not_tracking')
                continue

 
        input = {}
        if cal2 == 'on':
            if mission == 'jason3':
# Application of CAL-2 where known
                pathcal=os.path.join(os.path.abspath(os.path.dirname(__file__)), '../data/cal2/J3_MeanFilterKu')
                filter = np.loadtxt(pathcal)
                filter_norm = filter / np.mean(filter[11:115])
                input['waveform'] = S_waveform[
                    index_waveforms_row, index_waveforms_col, :] / filter_norm[11:115]
            elif mission == 'jason2':
                pathcal=os.path.join(os.path.abspath(os.path.dirname(__file__)), '../data/cal2/J2_MeanFilterKu')
                filter = np.loadtxt(pathcal)
                filter_norm = filter / np.mean(filter[11:115])
                input['waveform'] = S_waveform[
                    index_waveforms_row, index_waveforms_col, :] / filter_norm[11:115]
            elif mission == 'jason1':
                pathcal=os.path.join(os.path.abspath(os.path.dirname(__file__)), '../data/cal2/J1_MeanFilterKu')
                filter = np.loadtxt(pathcal)
                filter_norm = filter / np.mean(filter[11:115])
                input['waveform'] = S_waveform[
                    index_waveforms_row, index_waveforms_col, :] / filter_norm[11:115]
            elif mission == 'saral':
                pathcal=os.path.join(os.path.abspath(os.path.dirname(__file__)), '../data/cal2/ALK_MeanFilter')                
                filter = np.loadtxt(pathcal)
                filter_norm = filter / np.mean(filter)
                input['waveform'] = S_waveform[index_waveforms_row,
                                    index_waveforms_col, :] / filter_norm
            elif mission=='envisat':
                ' waveform '
                input['waveform'] = S_waveform[index_waveforms_row,:]  
            elif mission=='sentinel6_lrm':
                ' waveform '
                input['waveform'] = S_waveform[index_waveforms_row,:]                  
            elif mission=='cs2_lrm':
                ' waveform '
                input['waveform'] = S_waveform[index_waveforms_row,:]
            elif mission=='ers2':
                ' waveform '
                input['waveform'] = S_waveform[index_waveforms_row,index_waveforms_col,:]
            elif mission=='ers1':
                ' waveform '
                input['waveform'] = S_waveform[index_waveforms_row,index_waveforms_col,:]
            elif mission == 'envisat':
                input['waveform'] = S_waveform[index_waveforms_row, :]
            elif mission == 'cs2_lrm':
                input['waveform'] = S_waveform[index_waveforms_row, :]
            elif mission == 'swot': # uses J3 cal2 ... 
                pathcal=os.path.join(os.path.abspath(os.path.dirname(__file__)), '../data/cal2/J3_MeanFilterKu')
                filter = np.loadtxt(pathcal)
                filter_norm = filter / np.mean(filter[11:115])
                input['waveform'] = S_waveform[
                    index_waveforms_row, index_waveforms_col, :] / filter_norm[11:115]
#           
#                input['waveform'] = S_waveform[index_waveforms_row, index_waveforms_col, :] 
        else:
            input['waveform'] = S_waveform[index_waveforms_row,
                                index_waveforms_col, :]

        ' raw range in [m] '
        input['uralt'] = S_tracker[index_waveforms_row, index_waveforms_col]

        ' hsat '
        input['hsat'] = S_height[index_waveforms_row, index_waveforms_col]
        ' mission '
        input['mission'] = mission
        input['weights_type']  = weights_type
        input['costfunction']  = costfunction
        if mission=='ers2':
            input['mission'] = 'ers2_r_2cm'
        elif mission=='ers1':
            input['mission'] = 'ers1'    
        elif mission=='sentinel6_lrm':
            input['mission'] = 'sentinel6_lrm'
                    

        ' off nadir angle in degree '
        xifile=S_offnadir[index_waveforms_row, index_waveforms_col]
        xi=xifile        
        if (xifile > 32000): 
           xi=0.
        input['xi'] = xi
        if import_weights == 'yes':
            input['weights_flag'] = flag_edges
            input['weights'] = residual_std

        input['tau'] = tau
        input['Theta']  = Theta 
        input['SigmaP']  = SigmaP
        input['sqrtn']  = np.sqrt(nump)
        input['total_gate_number']=total_gate_number
        input['nominal_tracking_gate']=nominal_tracking_gate
        input['weight_outsub']  = weight_outsub
        input['smooth'] = smooth
#
# Calls retracker 
#
        retracker = WHALES_withRangeAndEpoch(input)
#
# Additional pass for SWH > 8 m ! activated with -S = ... 
# first we estimate a "recent wave height" from past waveforms
#
        SWH1D=0.
        if (np.shape(S_time)[1] >= 20):
            if (index_waveforms_row >0):
                SWH1D=np.nanmedian(swh_WHALES[index_waveforms_row-1,:])
        if (np.shape(S_time)[1] == 1):
            if (index_waveforms_row >=20):
                SWH1D=np.nanmedian(swh_WHALES[index_waveforms_row-20:index_waveforms_row,:])
# second, if this recent wave height exceeds a threshold, we retrack with smoothing as given by -S 
        if ((smooth_above_8m > 0) & (SWH1D > smooth_SWH_val ) & (SWH1D < 30. )):
                    input['smooth'] = smooth_above_8m
                    input['weights_type'] = 2
                    retracker = WHALES_withRangeAndEpoch(input)
#
# Post-processing
#
        # Quality flag of WHALES, based on the normalised fitting error on the leading edge
        if (retracker.Error) > 0.3 and (np.isnan(retracker.Error) == 0):
            Err_WHALES[index_waveforms_row, index_waveforms_col] = 1
        elif retracker.Error <= 0.3:
            Err_WHALES[index_waveforms_row, index_waveforms_col] = 0

        swh_WHALES[index_waveforms_row, index_waveforms_col] = retracker.SWH
#
# Added by Fabrice to keep track of actual part of waveform retracked 
#
        startgate_WHALES[index_waveforms_row, index_waveforms_col] = retracker.gate1
        endgate_WHALES[index_waveforms_row, index_waveforms_col] = retracker.gate2
        finalgate_WHALES[index_waveforms_row, index_waveforms_col] = retracker.gate3
        Epoch_WHALES[index_waveforms_row, index_waveforms_col] = retracker.Epoch
        noise_WHALES[index_waveforms_row, index_waveforms_col] = retracker.noise
        scale_WHALES[index_waveforms_row, index_waveforms_col] = retracker.scale

        swh_WHALES[index_waveforms_row,index_waveforms_col]               =retracker.SWH
        
        Epoch_WHALES[index_waveforms_row,index_waveforms_col]               =retracker.Epoch

# WARNING: NEED TO CHECK FOR CFOSAT AND SWOT ... 
        if mission in ['envisat','envisat_over']:
            sigma0_WHALES[index_waveforms_row,index_waveforms_col]            =retracker.Amplitude+S_atmos_corr[index_waveforms_row,index_waveforms_col]+ S_scaling_factor[index_waveforms_row,index_waveforms_col] -33.1133
        elif mission in ['sentinel6_lrm']:
            sigma0_WHALES[index_waveforms_row,index_waveforms_col]            =retracker.Amplitude+sig0_scaling_factor[index_waveforms_row,index_waveforms_col]+ 10 * np.log10(S_waveform_scale_factor[index_waveforms_row,index_waveforms_col])
        elif mission in ['ers2']: #The constant bias of 103.92 has been derived from the personal communication of David Brockley and it is also mentioned in the SGDR user manual for REAPER
            sigma0_WHALES[index_waveforms_row,index_waveforms_col]            =retracker.Amplitude+S_atmos_corr[index_waveforms_row,index_waveforms_col]+ S_scaling_factor[index_waveforms_row,index_waveforms_col] -103.92
        elif mission in ['ers1']: #The constant bias  has been derived from the personal communication of David Brockley and it is also mentioned in the SGDR user manual for REAPER
            sigma0_WHALES[index_waveforms_row,index_waveforms_col]            =retracker.Amplitude+S_atmos_corr[index_waveforms_row,index_waveforms_col]+ S_scaling_factor[index_waveforms_row,index_waveforms_col] -108.33    
        elif mission in ['jason2','jason1','saral','saral_igdr','jason3','swot']:
            sigma0_WHALES[index_waveforms_row,index_waveforms_col]            =retracker.Amplitude+S_atmos_corr[index_waveforms_row,index_waveforms_col]+ S_scaling_factor[index_waveforms_row,index_waveforms_col]
        elif mission in ['cs2_lrm']:
            sigma0_WHALES[index_waveforms_row,index_waveforms_col]            =retracker.Amplitude
        
        range_WHALES[index_waveforms_row,index_waveforms_col]             =retracker.range
        
        Amplitude_WHALES[index_waveforms_row,index_waveforms_col] = retracker.Norm_Amplitude

        range_WHALES[index_waveforms_row, index_waveforms_col] = retracker.range

        Amplitude_WHALES[
            index_waveforms_row, index_waveforms_col] = retracker.Norm_Amplitude

        # APPLICATION OF INSTRUMENTAL CORRECTION FOR SWH           
        if add_instr_corr_SWH == 'yes':
            if counting_swh == 0:
                swh_WHALES_instr_corr[
                    index_waveforms_row, index_waveforms_col], interpolator_instr_corr_SWH = compute_instr_corr_SWH_WHALES(
                    swh_WHALES[index_waveforms_row, index_waveforms_col],
                    my_path_instr_corr_SWH, mission)
                counting_swh = 1
            else:
                swh_WHALES_instr_corr[
                    index_waveforms_row, index_waveforms_col], interpolator_instr_corr_SWH = compute_instr_corr_SWH_WHALES(
                    swh_WHALES[index_waveforms_row, index_waveforms_col],
                    my_path_instr_corr_SWH, mission,
                    interpolator_instr_corr_SWH)
        if debug=='1':
            w_nc_fid.variables['normalized_waveform'][index_waveforms_row,index_waveforms_col,:] = retracker.D
            w_nc_fid.variables['weights'][index_waveforms_row,index_waveforms_col,:] = retracker.this_weights
            w_nc_fid.variables['fitted_waveform'][index_waveforms_row,index_waveforms_col,:] = retracker.model
    print('Hs:',np.shape(swh_WHALES),np.mean(swh_WHALES[index_waveforms_row,:]),' std:',np.std(swh_WHALES[index_waveforms_row,:]))


# now writes retracked values 
w_nc_var = w_nc_fid.createVariable('swh_WHALES_20hz', 'f4', ('time', 'records'),
                                   zlib=True)
w_nc_var.setncatts({'long_name': u"swh_WHALES_20hz", \
                    'units': u"m", \
                    'comment': u" "})
w_nc_fid.variables['swh_WHALES_20hz'][:] = swh_WHALES

w_nc_var = w_nc_fid.createVariable('swh_WHALES_instr_corr_20hz', 'f4',
                                   ('time', 'records'), zlib=True)
w_nc_var.setncatts({'long_name': u"swh_WHALES_instr_corr_20hz", \
                    'units': u"m", \
                    'comment': u" "})
w_nc_fid.variables['swh_WHALES_instr_corr_20hz'][:] = swh_WHALES_instr_corr

w_nc_var = w_nc_fid.createVariable('sigma0_WHALES_20hz', 'f4',
                                   ('time', 'records'), zlib=True)
w_nc_var.setncatts({'long_name': u"sigma0_WHALES_20hz", \
                    'units': u"dB", \
                    'comment': u" "})
w_nc_fid.variables['sigma0_WHALES_20hz'][:] = sigma0_WHALES

w_nc_var = w_nc_fid.createVariable('range_WHALES_20hz', 'f4',
                                   ('time', 'records'), zlib=True)
w_nc_var.setncatts({'long_name': u"range_WHALES_20hz", \
                    'units': u"m", \
                    'comment': u" "})
w_nc_fid.variables['range_WHALES_20hz'][:] = range_WHALES

w_nc_var = w_nc_fid.createVariable('epoch_WHALES_20hz', 'f4',
                                   ('time', 'records'), zlib=True)
w_nc_var.setncatts({'long_name': u"epoch_WHALES_20hz", \
                    'units': u"m", \
                    'comment': u" "})
w_nc_fid.variables['epoch_WHALES_20hz'][:] = Epoch_WHALES

w_nc_var = w_nc_fid.createVariable('swh_WHALES_qual_20hz', 'i2',
                                   ('time', 'records'), zlib=True)
w_nc_var.setncatts({'long_name': u"quality flag for Significant waveheight", \
                    'units': u"count", \
                    'comment': u"0=Good, 1=Bad"})
w_nc_fid.variables['swh_WHALES_qual_20hz'][:] = Err_WHALES

#
# Added by Fabrice to keep track of range gate indices that define the subwaveform 
#
w_nc_var = w_nc_fid.createVariable('startgate_WHALES', 'i2',
                                   ('time', 'records'), zlib=True)
w_nc_var.setncatts({'long_name': u"first index of subwaveform"})
w_nc_fid.variables['startgate_WHALES'][:] = startgate_WHALES

w_nc_var = w_nc_fid.createVariable('endgate_WHALES', 'i2',
                                   ('time', 'records'), zlib=True)
w_nc_var.setncatts({'long_name': u"last index of subwaveform"})
w_nc_fid.variables['endgate_WHALES'][:] = endgate_WHALES
w_nc_var = w_nc_fid.createVariable('finalgate_WHALES', 'i2',
                                   ('time', 'records'), zlib=True)
w_nc_var.setncatts({'long_name': u"last index of retracking"})
w_nc_fid.variables['finalgate_WHALES'][:] = finalgate_WHALES

w_nc_var = w_nc_fid.createVariable('scale_WHALES', 'f4',
                                   ('time', 'records'), zlib=True)
w_nc_var.setncatts({'long_name': u"scaling for waveform normalization"})
w_nc_fid.variables['scale_WHALES'][:] = scale_WHALES
w_nc_var = w_nc_fid.createVariable('noise_WHALES', 'f4',
                                   ('time', 'records'), zlib=True)
w_nc_var.setncatts({'long_name': u"noise level of waveform"})
w_nc_fid.variables['noise_WHALES'][:] = noise_WHALES



w_nc_fid.close()  # close the new file
