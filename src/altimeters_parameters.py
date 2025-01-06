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
import netCDF4
from netCDF4 import Dataset


######################  Defines parameters
def  instrument_parameters(mission)  :
# Default values for the Jason series 
    tau=3.125 #gate spacing in ns
    SigmaP=0.513*tau
    nump=90
    total_gate_number=104
    nominal_tracking_gate=31     
                
# Mission-dependent parameters and files to be loaded
    if mission.lower() in ['ers2','ers2_r','ers2_r_2cm']:           
        tau=3.03
        Theta=1.3 *np.pi/180
        SigmaP=0.513*tau    
        nump=50
        total_gate_number=64
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
        my_path_weights = '../data/weights/weights_n1.mat'
    elif mission in ['ers2','ers1'] :
        my_path_instr_corr_SWH='instr_corr/SWHinstrcorr_WHALES_jason3SGDRd.mat' 
        my_path_weights='../data/weights/weights_ers2.pkl'       
    elif mission in ['jason1']:
        my_path_instr_corr_SWH = 'instr_corr/SWHinstrcorr_MLE4_jason1SGDRc.mat'
        my_path_weights = '../data/weights/weights.mat'
    elif mission in ['jason2']:
        my_path_instr_corr_SWH = 'instr_corr/SWHinstrcorr_WHALES_jason2SGDRd.mat'
        my_path_weights = '../data/weights/weights_J2.pkl'
    elif mission in ['jason3', 'jason3f','jason3f2','swot']:
        my_path_instr_corr_SWH = 'instr_corr/SWHinstrcorr_WHALES_jason3SGDRd.mat'
        my_path_weights = '../data/weights/weights_J2.pkl'
    elif mission.lower() in ['altika', 'saral', 'saral_igdr']:
        my_path_instr_corr_SWH = ''
        my_path_weights = '../data/weights/weights_alt.mat'
    elif mission in ['cs2_lrm']:
        my_path_instr_corr_SWH = ''
        my_path_weights = '../data/weights/weights_cs2_lrm.mat'
    return my_path_instr_corr_SWH,my_path_weights
    
    ######################  Generic reading of altimeter data: retracked parameters
def  alti_read_l2lr(mission,filename):
    '''
    reads altimeter data (LRM 1Hz only) from file name. 
    The outout is a xarray dataset 
    '''

    import netCDF4
    import xarray as xr
    from netCDF4 import Dataset  
    S = netCDF4.Dataset(filename, 'r')
    if mission.lower() in ['jason1','jason2','jason3']:
    # example file='JA2_GPS_2PdP011_200_20081026_233206_20081027_002819.nc'
        swh1 = np.ma.getdata(S.variables['swh_ku'][:])   # this is MLE4
        lat1  = np.ma.getdata(S.variables['lat'][:])
        lon1  = np.ma.getdata(S.variables['lon'][:])
        time1 = np.ma.getdata(S.variables['time'][:])
        flag1 = np.ma.getdata(S.variables['qual_alt_1hz_swh_ku'][:])
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

    if mission.lower() in ['jason3f']:
    # example file='JA2_GPS_2PdP011_200_20081026_233206_20081027_002819.nc'
        swh1 = np.ma.getdata(S['data_01']['ku'].variables['swh_ocean'][:])   # this is MLE4
        lat1  = np.ma.getdata(S['data_01'].variables['lat'][:])
        lon1  = np.ma.getdata(S['data_01'].variables['lon'][:])
        time1 = np.ma.getdata(S['data_01'].variables['time'][:])
        flag1 = np.ma.getdata(S['data_01']['ku'].variables['swh_ocean_mle3_compression_qual'][:])
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 


    if mission in ['saral', 'altika']:
    # example file='SRL_GPS_2PfP001_0641_20130405_141055_20130405_150113.CNES.nc'
        swh1 = np.ma.getdata(S.variables['swh'][:])   # this is MLE4?
        lat1  = np.ma.getdata(S.variables['lat'][:])
        lon1  = np.ma.getdata(S.variables['lon'][:])
        time1 = np.ma.getdata(S.variables['time'][:])
        flag1 = np.ma.getdata(S.variables['qual_alt_1hz_swh'][:])
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

    if mission.lower() in ['cfosat']:
    # example file='CFO_OP07_SWI_L2_____F_20241221T043021_20241221T054444.nc'
        S = netCDF4.Dataset(filename, 'r')
        swh1 = np.ma.getdata(S.variables['nadir_swh_1Hz'][:])
    #    swh_sgdr2 = np.ma.getdata(S.variables['nadir_swh_native'][:])
    #S_waveform = np.ma.getdata(S.variables['waveforms_40hz'][:])
        lat1  = np.ma.getdata(S.variables['lat_nadir_1Hz'][:])
        lon1  = np.ma.getdata(S.variables['lon_nadir_1Hz'][:])
#        time1 = np.ma.getdata(S.variables['time_nadir_1Hz'][:])
        time1 = S.variables['time_nadir_1Hz'][:]
        flag1 = np.ma.getdata(S.variables['flag_valid_swh_1Hz'][:])
        timeref= "2009-01-01 00:00:00 0:00"
    #reference_time = pd.Timestamp("2014-09-05")
    ds = xr.Dataset(
        {   "swh_1hz": (["time"], swh1),
            "lon_1hz": (["time"], lon1),
            "lat_1hz": (["time"], lat1),
            "flag_1hz": (["time"], flag1),
        },
        coords={
            "time": time1,
            "reference_time": timeref,
        },
        )
    return ds


######################  Generic reading of altimeter data: retracked 1Hz parameters
def  alti_read_l2lr(mission,filename):
    '''
    reads altimeter data (LRM 1Hz only) from file name. 
    The outout is a xarray dataset 
    '''

    import xarray as xr

    S = netCDF4.Dataset(filename, 'r')
    if mission.lower() in ['ers2','ers2_r_2cm']:
    # example file='CS_LTA__SIR_LRM_2__20101231T000154_20101231T000927_E001.nc'
        swh1 = np.ma.getdata(S.variables['swh'][:])   # this is MLE4
        lat1  = np.ma.getdata(S.variables['lat'][:])
        lon1  = np.ma.getdata(S.variables['lon'][:])
        time1 = np.ma.getdata(S.variables['time'][:])
        try:
            flag1 = 3 - np.ma.getdata(S.variables['swh_quality_level'][:])
        except KeyError:
           print("Variable 'swh_quality_level' not found in the dataset.")
           # Handle missing data or use a default value
           flag1 = np.zeros_like(S.variables['time'][:])  # or other default behavior

        timeref= "1981-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 
    if mission.lower() in ['cryosat2']:
    # example file='CS_LTA__SIR_LRM_2__20101231T000154_20101231T000927_E001.nc'
        swh1 = np.ma.getdata(S.variables['swh_ku'][:])   # this is MLE4
        lat1  = np.ma.getdata(S.variables['lat'][:])
        lon1  = np.ma.getdata(S.variables['lon'][:])
        time1 = np.ma.getdata(S.variables['time'][:])
        flag1 = np.ma.getdata(S.variables['qual_alt_1hz_swh_ku'][:])
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

    if mission.lower() in ['jason1','jason2','jason3']:
    # example file='JA2_GPS_2PdP011_200_20081026_233206_20081027_002819.nc'
        swh1 = np.ma.getdata(S.variables['swh_ku'][:])   # this is MLE4
        lat1  = np.ma.getdata(S.variables['lat'][:])
        lon1  = np.ma.getdata(S.variables['lon'][:])
        time1 = np.ma.getdata(S.variables['time'][:])
        flag1 = np.ma.getdata(S.variables['qual_alt_1hz_swh_ku'][:])
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

    if mission.lower() in ['jason3f','swot']:
    # example file='JA2_GPS_2PdP011_200_20081026_233206_20081027_002819.nc'
        swh1 = np.ma.getdata(S['data_01']['ku'].variables['swh_ocean'][:])   # this is MLE4
        lat1  = np.ma.getdata(S['data_01'].variables['latitude'][:])
        lon1  = np.ma.getdata(S['data_01'].variables['longitude'][:])
        time1 = np.ma.getdata(S['data_01'].variables['time'][:])
        flag1 = np.ma.getdata(S['data_01']['ku'].variables['swh_ocean_compression_qual'][:])
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

    if mission in ['saral', 'altika']:
    # example file='SRL_GPS_2PfP001_0641_20130405_141055_20130405_150113.CNES.nc'
        swh1 = np.ma.getdata(S.variables['swh'][:])   # this is MLE4?
        lat1  = np.ma.getdata(S.variables['lat'][:])
        lon1  = np.ma.getdata(S.variables['lon'][:])
        time1 = np.ma.getdata(S.variables['time'][:])
        flag1 = np.ma.getdata(S.variables['qual_alt_1hz_swh'][:])
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

    if mission.lower() in ['cfosat']:
    # example file='CFO_OP07_SWI_L2_____F_20241221T043021_20241221T054444.nc'
        S = netCDF4.Dataset(filename, 'r')
        swh1 = np.ma.getdata(S.variables['nadir_swh_1Hz'][:])
    #    swh_sgdr2 = np.ma.getdata(S.variables['nadir_swh_native'][:])
    #S_waveform = np.ma.getdata(S.variables['waveforms_40hz'][:])
        lat1  = np.ma.getdata(S.variables['lat_nadir_1Hz'][:])
        lon1  = np.ma.getdata(S.variables['lon_nadir_1Hz'][:])
#        time1 = np.ma.getdata(S.variables['time_nadir_1Hz'][:])
        time1 = S.variables['time_nadir_1Hz'][:]
        flag1 = np.ma.getdata(S.variables['flag_valid_swh_1Hz'][:])
        timeref= "2009-01-01 00:00:00 0:00"
    #reference_time = pd.Timestamp("2014-09-05")
    ds = xr.Dataset(
        {   "swh_1hz": (["time"], swh1),
            "lon_1hz": (["time"], lon1),
            "lat_1hz": (["time"], lat1),
            "flag_1hz": (["time"], flag1),
        },
        coords={
            "time": time1,
            "reference_time": timeref,
        },
        )
    return ds
    
######################  Generic reading of altimeter data: retracked 1Hz parameters
def  alti_read_l2lr_cci(mission,filename):
    '''
    reads altimeter data (LRM 1Hz only) from file name. 
    The outout is a xarray dataset 
    '''

    import xarray as xr
    S = netCDF4.Dataset(filename, 'r')
    swh1 = np.ma.getdata(S.variables['swh'][:])   
    rms1  = np.ma.getdata(S.variables['swh_rms'][:])
    lat1  = np.ma.getdata(S.variables['lat'][:])
    lon1  = np.ma.getdata(S.variables['lon'][:])
    time1 = np.ma.getdata(S.variables['time'][:])
    flag1 = 3-np.ma.getdata(S.variables['swh_quality_level'][:])
    timeref= "1981-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 


    ds = xr.Dataset(
        {   "swh_1hz": (["time"], swh1),
            "swh_rms": (["time"], rms1),
            "lon_1hz": (["time"], lon1),
            "lat_1hz": (["time"], lat1),
            "flag_1hz": (["time"], flag1),
        },
        coords={
            "time": time1,
            "reference_time": timeref,
        },
        )
    return ds

######################  Generic reading of altimeter data: waveforms and 20Hz parameters
def  alti_read_l2hrw(mission,filename):
    '''
    reads altimeter data (LRM 1Hz only) from file name. 
    The outout is a xarray dataset 
    '''
    import xarray as xr
    S = netCDF4.Dataset(filename, 'r')
    if mission.lower() in ['ers1','ers2']:
        #Time at 1-Hz for interpolation of fields available only at 1-Hz
        time1=np.ma.getdata( S.variables['time'][:] )
        S_time = np.ma.getdata(S.variables['time_20hz'][:]).flatten()
        nhf=20
        n1=len(time1)
        nall=len(S_time)
        nlr=nall//nhf
        nal=nlr*nhf

        S_time=np.ma.getdata( S.variables['time_20hz'][:] )
        #S_time=np.reshape(S_time,(np.shape(S_time)[0],1) )
    
    
        height1=np.ma.getdata( S.variables['alt_20hz'][:] )
        swh1=np.ma.getdata( S.variables['swh_20hz'][:] )
        tracker1=np.ma.getdata( S.variables['tracker_range_20hz'][:] )
        range1=np.ma.getdata( S.variables['ocean_range_20hz'][:] )
        waveforms=np.ma.getdata( S.variables['ku_wf'][:] ).astype(np.float64)

        lat1=np.ma.getdata( S.variables['lat_20hz'][:] )
        lon1=np.ma.getdata( S.variables['lon_20hz'][:] )
        off=np.ma.getdata( S.variables['off_nadir_angle_wf_20hz'][:] )
        off=np.zeros_like(off) # THE OFF NADIR FIELD IN ERS2 IS EMPTY!!!!!

        atmos_corr=np.ma.getdata( S.variables['atmos_corr_sig0'][:] )
        atmos_corr=np.transpose(np.tile(atmos_corr,(np.shape(S_time)[1],1)))
        scaling_factor=np.ma.getdata( S.variables['scaling_factor_20hz'][:] )
        timeref= "1981-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

        
    if mission.lower() in ['jason1','jason2','jason3']:
    # example file='JA2_GPS_2PdP011_200_20081026_233206_20081027_002819.nc'
        time1 = np.ma.getdata(S.variables['time'][:])
        S_time = np.ma.getdata(S.variables['time_20hz'][:]).flatten()
        nhf=20
        n1=len(time1)
        nall=len(S_time)
        nlr=nall//nhf
        nal=nlr*nhf
        swh1 = np.ma.getdata(S.variables['swh_20hz_ku'][:])  # this is MLE4
        lat1= np.ma.getdata(S.variables['lat_20hz'][:])
#,np.shape(S_time),np.shape(S_swh),np.shape(S_lat))
        #swh1 = np.reshape(S_swh[0:nal], (nlr,nhf))
        #lat1 = np.reshape(S_lat[0:nal], (nlr,nhf))
        lon1= np.ma.getdata(S.variables['lon_20hz'][:])
        #lon1 =np.reshape(S_lon[0:nal], (nlr,nhf))

                
        waveforms = np.ma.getdata(S.variables['waveforms_20hz_ku'][:])
        #waveforms = np.reshape(
        #S_waveform[0:nal,:], (nlr, nhf, np.shape(S_waveform)[1]))

        flag1 = np.ma.getdata(S.variables['qual_alt_1hz_swh_ku'][:])
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

    if mission.lower() in ['jason3f','swot']:
        # AVISO SGDR version F
        time1 = np.ma.getdata(S['data_01'].variables['time'][:])
        S_time = np.ma.getdata(S['data_20'].variables['time'][:])
        nhf=20
        n1=len(time1)
        nall=len(S_time)
        nlr=nall//nhf
        nal=nlr*nhf
        
        #print('time:',time1[0],'##',S_time[0:20]-time1[0])
        #print('time:',time1[nlr-1],'##',S_time[(nlr-1)*nhf:nlr*nhf]-time1[nlr-1])
        S_time = np.reshape(S_time[0:nal], (nlr,nhf))   # warning this reshaping throws away the last few 20 Hz sample 
        #print('time shape:',n1,S_time.shape)
       
        S_swh = np.ma.getdata(S['data_20']['ku'].variables['swh_ocean'][:]) # this is MLE4
        swh1 = np.reshape(S_swh[0:nal], (nlr,nhf))

        S_lat= np.ma.getdata(S['data_20'].variables['latitude'][:])
        lat1 = np.reshape(S_lat[0:nal], (nlr,nhf))

        S_lon= np.ma.getdata(S['data_20'].variables['longitude'][:])
        lon1 =np.reshape(S_lon[0:nal], (nlr,nhf))
        
        
        S_waveform = np.ma.getdata(
        S['data_20']['ku'].variables['power_waveform'][:])
        waveforms = np.reshape(
        S_waveform[0:nal,:], (nlr, nhf, np.shape(S_waveform)[1]))
        
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

    if mission in ['saral', 'altika']:
    # example file='SRL_GPS_2PfP001_0641_20130405_141055_20130405_150113.CNES.nc'
        swh2 = np.ma.getdata(S.variables['swh_40hz'][:])   # this is MLE4?
        lat2  = np.ma.getdata(S.variables['lat_40hz'][:])
        lon2  = np.ma.getdata(S.variables['lon_40hz'][:])
        time2 = np.ma.getdata(S.variables['time_40hz'][:])
        time1 = np.ma.getdata(S.variables['time'][:])
        off2 = np.ma.getdata(S.variables['off_nadir_angle_pf'][:])
        waveforms = np.ma.getdata(S.variables['waveforms_40hz'][:])

        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

    if mission.lower() in ['cfosat']:
    # example file='CFO_OP07_SWI_L2_____F_20241221T043021_20241221T054444.nc'
        S = netCDF4.Dataset(filename, 'r')
        swh1 = np.ma.getdata(S.variables['nadir_swh_1Hz'][:])
    #    swh_sgdr2 = np.ma.getdata(S.variables['nadir_swh_native'][:])
        waveforms = np.ma.getdata(S.variables['waveforms_40hz'][:])
        lat1  = np.ma.getdata(S.variables['lat_nadir_1Hz'][:])
        lon1  = np.ma.getdata(S.variables['lon_nadir_1Hz'][:])
#        time1 = np.ma.getdata(S.variables['time_nadir_1Hz'][:])
        time1 = S.variables['time_nadir_1Hz'][:]
        flag1 = np.ma.getdata(S.variables['flag_valid_swh_1Hz'][:])
        timeref= "2009-01-01 00:00:00 0:00"
    #reference_time = pd.Timestamp("2014-09-05")
    ds = xr.Dataset(
        {   "swh2d": (["time","meas_ind"], swh1),
            "lon2d": (["time","meas_ind"], lon1),
            "lat2d": (["time","meas_ind"], lat1),
            "waveforms": (["time","meas_ind","wvf_ind"], waveforms),
        },
        coords={
            "time": time1[0:nlr],
            "reference_time": timeref,
        },
        )
    return ds
    
    ######################  Generic reading of altimeter data: waveforms and 20Hz parameters
def  alti_read_l2hrw_cci(mission,filename):
    '''
    reads altimeter data (LRM 1Hz only) from file name. 
    The outout is a xarray dataset 
    '''

    S = netCDF4.Dataset(filename, 'r')
    if mission.lower() in ['jason1','jason2','jason3']:
    # example file='JA2_GPS_2PdP011_200_20081026_233206_20081027_002819.nc'
        swh1 = np.ma.getdata(S.variables['swh_ku'][:])   # this is MLE4
        lat1  = np.ma.getdata(S.variables['lat'][:])
        lon1  = np.ma.getdata(S.variables['lon'][:])
        time1 = np.ma.getdata(S.variables['time'][:])
        flag1 = np.ma.getdata(S.variables['qual_alt_1hz_swh_ku'][:])
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 

    if mission.lower() in ['jason3f']:
        # AVISO SGDR version F
        S_time = np.ma.getdata(S.variables['time_20hz'][:])
        nhf=20
        nall=len(S_time)
        nlr=nall//nhf
        nal=nlr*nhf
        
        #print('time:',time1[0],'##',S_time[0:20]-time1[0])
        #print('time:',time1[nlr-1],'##',S_time[(nlr-1)*nhf:nlr*nhf]-time1[nlr-1])
        S_time = np.reshape(S_time[0:nal], (nlr,nhf))   # warning this reshaping throws away the last few 20 Hz sample 
        time1 = np.mean(S_time,axis=1)
        
        S_swh = np.ma.getdata(S.variables['swh_WHALES_20hz'][:]) # this is MLE4
        swh1 = np.reshape(S_swh[0:nal], (nlr,nhf))

        S_lat= np.ma.getdata(S.variables['lat_20hz'][:])
        lat1 = np.reshape(S_lat[0:nal], (nlr,nhf))

        S_lon= np.ma.getdata(S.variables['lon_20hz'][:])
        lon1 =np.reshape(S_lon[0:nal], (nlr,nhf))
        
        
#        S_waveform = np.ma.getdata(
#        S['data_20']['ku'].variables['power_waveform'][:])
#        waveforms = np.reshape(
#        S_waveform[0:nal,:], (nlr, nhf, np.shape(S_waveform)[1]))
        
        timeref= "2000-01-01 00:00:00.0"			# WARNING: this should be read from the attribute of the time variable ... 


    ds = xr.Dataset(
        {   "swh2d": (["time","meas_ind"], swh1),
            "lon2d": (["time","meas_ind"], lon1),
            "lat2d": (["time","meas_ind"], lat1),
 #           "waveforms": (["time","meas_ind","wvf_ind"], waveforms),
        },
        coords={
            "time": time1[0:nlr],
            "reference_time": timeref,
        },
        )
    return ds


################################################################################################################    
def  processing_choices(mission)  :
    thra=0.1   # threshold for normalized waveform at second range gate of leading edge
    # This second threshold was introduced by FA (to recover previous behavior, set thrb to 0) . 
    thrb=0 #0.7   # threshold for lowest normalized waveform value beyond which the leading edge may stop. 
    minHs=0.2
    maxHs=30
    noisemin=0 #1
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

   
    elif mission.lower() == 'jason2' or mission.lower() == 'jason1' or mission.lower() == 'jason3' or mission.lower() == 'jason3f' or mission.lower() == 'jason3f2' or mission.lower() == 'swot': 
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
                
    elif mission.lower() in ['ers2','ers2_r_2cm']:
                index_originalbins=np.arange(0,63,1) 
                total_gate_number=64
                noisegates=np.arange(10,13); #gates used to estimate Thermal Noise (see mail from David Brockley)
                tau=3.03 #gate width in nanoseconds
                startgate=8 #First gate to be considered in the retracking window, we choose this because first gates are often corrupted
                ALEScoeff0=8.90 #experimental values for SWH. it is the constant term in the definition of the number of gates to be considered in the retracking
                                #after the middle of the leading edge
                ALEScoeff1=2.03 #This is the slope of the WHALES relationship between tolerance of precision and width of the subwaveform   
                Err_tolerance_vector=0.3; #Tolerance on the (normalised) fitting error of the waveform. It can be used, for example,
                                                        #to retrack the same waveform in a different way if fitting performances are not satisfactory  
                # # AKIMA INTERPOLATION                
                # waveform_resampled=np.empty(np.size(waveform))*np.nan
                # y=waveform # 11 May 2015: forced conversion to double                
                # x=np.arange(1*tau,64*tau+0.01,tau)
                # x[-1]=round(x[-1],5)
                # xi=np.arange(1*tau,64*tau+0.01,tau/2)
                # xi[-1]=round(xi[-1],5)
                # yi=self.akima_interpolate(x,y,xi)
                # waveform_resampled=np.append(yi,yi[-1]) # repeat last sample to make 208                
                # # END AKIMA INTERPOLATION
                # waveform = waveform_resampled    
    else:
                print("unknown mission")
                sys.exit(0)
    return(noisegates,startgate,ALEScoeff0,ALEScoeff1,Err_tolerance_vector,thra,thrb,minHs,maxHs,noisemin) 

