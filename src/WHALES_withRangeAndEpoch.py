# -*- coding: utf-8 -*-
"""
Created on 11 February 2019

@author: passaro

WHALES Retracker corresponding to the WHALES version of the 4.12.2018, but with Epoch and Range in Output

The WHALES algorithm is described in 
Passaro, M., & Algorithm Development Team. (2021). Algorithm theoretical basis document (atbd), 
sea state climate change initiative (Tech. Rep.). European Space Agency. 
Retrieved from https://climate.esa.int/media/documents/Sea_State_cci_ATBD_v3.0-signed.pdf

Modification history:
	2023-01-10: corrections for SARAL by Marine De Carlo 
        2024-06-30: adaptation for wavesALTI package by F. Ardhuin 
        2024-07-19: added possibility to use 1/waveform for the weights
        2024-12-19: export of gate1, gate2 and gate3 indices, clean-up


"""


from Retracker_MP    import *
from altimeters_parameters import processing_choices
from waveform_models import waveform_brown_LS,waveform_brown_ML,wf_brown_eval
import scipy
from scipy           import stats
import cmath #handling of complex square root
from scipy           import optimize
from scipy.optimize  import minimize
from scipy           import special
from scipy           import signal
from math            import erf
import math
import numpy as np
import matplotlib.pyplot as plt
' WHALES Retracker '
import time



class WHALES_withRangeAndEpoch(Retracker_MP):
    
    def __init__(self,config):
        Retracker_MP.__init__(self, config)
        self.retrack_MP()


    def NM_fit(self,xdata,ydata,Zeta,tau,Theta,SigmaP,altitude,initial_conditions,mission,weights,weightflag,modelcost) :
        #Nelder-Mead fit of a waveform 
    
        #IT NEEDS:
            #1) xdata in NANOSECONDS (ex. xdata=[0*tau:1:127*tau])
            #2) ydata (waveform coefficient) in normalized units
            #3) Off-nadir angle in radians
    
        #OUTPUT:
            #1) x: vector of the estimated parameters
            #2) Wt: fitted waveform
            #3) exitflag: exit flag to check convergence
            #4) Err: Fitting Error on the portion of waveform given to the function
            #5) SWH: Significant Wave Height
        
        # WARNING: ZETA is xi converted to radians 
        
        incognita=initial_conditions
        
        c=3.0*(10**8) #Light speed
        H=altitude
        Ri=6378.1363*(10**3) #Earth radius

        Gamma=0.5 * (1/math.log(2))*np.sin(Theta)*np.sin(Theta) # antenna beamwidth parameter
     
        b_xi = np.cos (2*Zeta) - ((np.sin(2*Zeta))**2)/Gamma
        a=( (4/Gamma)*(c/H) * 1/(1+H/Ri))
        c_xi=b_xi* ( (4/Gamma)*(c/H) * 1/(1+H/Ri))
    
        a=a/1000000000 #/ns
        c_xi=c_xi/1000000000 #1/ns


        if modelcost == 'brown_LS':
           xopt = minimize(waveform_brown_LS, incognita, args=((ydata,Gamma,Zeta,xdata,SigmaP,c_xi,weights,weightflag),) ,method='Nelder-Mead',options={'disp': False})
        elif modelcost == 'brown_ML':
           xopt = minimize(waveform_brown_ML, incognita, args=((ydata,Gamma,Zeta,xdata,SigmaP,c_xi,weights,weightflag),) ,method='Nelder-Mead',options={'disp': False})

        
        if xopt.success == True:
            exitflag=1 #xopt[4]
        else:
            exitflag=0
        
        #Calculation of the fitted waveform
        x=xopt.x #These are the parameters estimated in minimisation

        SWH_squared=( - SigmaP**2 + x[1]**2 )
        if SWH_squared>=0 :    
            SWH=cmath.sqrt( np.abs(- SigmaP**2 + x[1]**2)  ).real * (2*0.3)   # factor 2*0.3 converts sigma_time in nanosecond to 4*sigma in meters
        else:
            SWH=- cmath.sqrt( np.abs(- SigmaP**2 + x[1]**2) ).real * (2*0.3)
        SigmaS=SWH/(2*0.3) 
        Sigma=cmath.sqrt(SigmaS**2+SigmaP**2).real
        t=xdata-x[0]
        v=c_xi*(t-c_xi*(Sigma**2)/2)
        u=np.divide( (t-c_xi*(Sigma**2)) , (np.sqrt(2)*Sigma) )
        A=x[2]/2*np.exp((-4/Gamma)*(np.sin(Zeta))**2)
        Wt=A*np.exp(-v)*(1+scipy.special.erf(u))
        
        if np.size(ydata)>0: 
            Err=np.sqrt( 1./np.size(ydata) * np.sum( (ydata-Wt)**2 ))  #Fitting Error
        else:
            Err=np.nan
        
        return x, Wt, exitflag, Err, SWH
    
    
    def Conversion_NMbrown(self,x,tau,nominal_tracking_gate) :
        
        # This function converts the Epoch estimated by NM_fit into an Epoch referred to the nominal point of the mission
        # In input it takes x, the vector of parameters estimated by NM_fit            
        c=0.3 # Light speed divided by a factor 10^-9
       
        t0=x[0] # ns
        Sigma=x[1]
        Au=x[2]
        #============== Calculate Epoch from t0=========================
        Epoch=t0 - self.nominal_tracking_gate*self.tau;  
        Epoch=Epoch*(c/2)  #m conversion from ns to meters(1gate=0.468m->c*tau/2)
        
        return(Epoch,t0,Sigma,Au)



    def akima_interpolate(self, x, y, x_new, axis=-1, out=None):
        """
        Return interpolated data using Akima's method.

        This Python implementation is inspired by the Matlab(r) code by
        N. Shamsundar. It lacks certain capabilities of the C implementation
        such as the output array argument and interpolation along an axis of a
        multidimensional data array.

        Parameters
        ----------
        x : array like
            1D array of monotonically increasing real values.
        y : array like
            N-D array of real values. y's length along the interpolation
            axis must be equal to the length of x.
        x_new : array like
            New independent variables.
        axis : int
            Specifies axis of y along which to interpolate. Interpolation
            defaults to the last axis of y.
        out : array
            Optional array to receive results. Dimension at the specified
            axis must equal the length of x.
        """
        x = np.array(x, dtype=np.float64, copy=True)
        y = np.array(y, dtype=np.float64, copy=True)
        xi = np.array(x_new, dtype=np.float64, copy=True)

        if axis != -1 or out is not None or y.ndim != 1:
            raise NotImplementedError("implemented in C extension module")

        if x.ndim != 1 or xi.ndim != 1:
            raise ValueError("x-arrays must be one dimensional")

        n = len(x)
        if n < 3:
            raise ValueError("array too small")
        if n != y.shape[axis]:
            raise ValueError("size of x-array must match data shape")

        dx = np.diff(x)
        if any(dx <= 0.0):
            raise ValueError("x-axis not valid")

        if any(xi < x[0]) or any(xi > x[-1]):
            raise ValueError("interpolation x-axis out of bounds")

        m = np.diff(y) / dx
        mm = 2.0 * m[0] - m[1]
        mmm = 2.0 * mm - m[0]
        mp = 2.0 * m[n - 2] - m[n - 3]
        mpp = 2.0 * mp - m[n - 2]

        m1 = np.concatenate(([mmm], [mm], m, [mp], [mpp]))

        dm = np.abs(np.diff(m1))
        f1 = dm[2:n + 2]
        f2 = dm[0:n]
        f12 = f1 + f2

        ids = np.nonzero(f12 > 1e-9 * np.max(f12))[0]
        b = m1[1:n + 1]

        b[ids] = (f1[ids] * m1[ids + 1] + f2[ids] * m1[ids + 2]) / f12[ids]
        c = (3.0 * m - 2.0 * b[0:n - 1] - b[1:n]) / dx
        d = (b[0:n - 1] + b[1:n] - 2.0 * m) / dx ** 2

        bins = np.digitize(xi, x)
        bins = np.minimum(bins, n - 1) - 1
        bb = bins[0:len(xi)]
        wj = xi - x[bb]

        return ((wj * d[bb] + c[bb]) * wj + b[bb]) * wj + y[bb]


 ###############################################################################################   
    def retrack_MP(self):

        mission = self.mission
        tau     = self.tau
        waveform = np.array(self.waveform)
        
            
        # IN THE FOLLOWING LINES; THE SPECIFIC CHARACTERISTICS OF EACH MISSIONS ARE DEFINED
        # NOTE THAT:
        # Waveforms are not oversampled, because Jason has been tested with the addition of 
        # weights, whose distribution might change if we oversample the waveform.   
        noisegates,startgate,ALEScoeff0,ALEScoeff1,Err_tolerance_vector,thra,thrb,minHs,maxHs,noisemin,interpolation_factor,index_originalbins,total_gate_number=processing_choices(mission)
        
        # INITIAL DEFINITION OF EXIT VARIABLES

        self.Epoch=np.nan
        self.SWH=np.nan
        self.SWH_yang=np.nan
        self.Amplitude=np.nan
        self.Error = np.nan                  
        self.range = np.nan                     
        self.uncssh = np.nan #Uncorrected Sea Surface Height: not computed 
        Epoch_LESfive=np.nan                
        if (interpolation_factor > 1): 
            # AKIMA INTERPOLATION                
            waveform_resampled=np.empty(np.size(waveform))*np.nan
            y=waveform # 11 May 2015: forced conversion to double           
            x=np.arange(1*tau,int(total_gate_number/interpolation_factor)*tau+0.01,tau)
            x[-1]=round(x[-1],5)
            xi=np.arange(1*tau,int(total_gate_number/interpolation_factor)*tau+0.01,tau/2)
            xi[-1]=round(xi[-1],5)
            yi=self.akima_interpolate(x,y,xi)
            waveform_resampled=np.append(yi,yi[-1]) # repeat last sample to make 208                
            # END AKIMA INTERPOLATION
            waveform = waveform_resampled                
                
                
# STEP 1: NOISE ESTIMATION                 
        estimated_noise = np.nanmean(waveform[noisegates])
        if noisemin==1:
            estimated_noise2=np.nanmin(waveform[noisegates[0]:self.nominal_tracking_gate])*0.99999
            if (estimated_noise > estimated_noise2*3):  # keeps the mean value if possible 
                estimated_noise=estimated_noise2
                noisegates[0]=np.argmin(waveform[noisegates[0]+4:self.nominal_tracking_gate])
        C = waveform - estimated_noise;

        self.noise=estimated_noise
# STEP 2: WAVEFORM NORMALISATION
        igoodsample=C>np.max([5, 2*estimated_noise]) #5 and 2 are arbitrary factors here   # 5 should be increased for ERS ... also need to exclue the noisegates ... 
        normalize_factor=1.3*np.nanmedian(C[igoodsample]);        
        D=C/normalize_factor
        self.scale=normalize_factor
            
# STEP 3: DEFINING RETRACKING ABSCISSA (xdata)
        xdata=np.empty([total_gate_number,])*np.nan
        xdata[index_originalbins]=np.arange(0,((np.size(index_originalbins)-1)*tau)+tau,tau)
        #xdata[index_originalbins]=np.arange(0,((np.size(index_originalbins))*tau),tau)
        iii=1            
        for iii in np.arange(1,total_gate_number,2) :
            if iii==total_gate_number-1 :
                xdata[iii]=xdata[iii-1]+(xdata[iii-1]-xdata[iii-2])
            else:
                xdata[iii]=(xdata[iii+1]+xdata[iii-1])/2
                  
            
        
# STEP 4: DETECTION OF LEADING EDGE
        self.Wt_all_yang= np.empty(total_gate_number)*np.nan       
        self.Wt_all_LESfive= np.empty(total_gate_number)*np.nan  
        
            # STEP 4: FIND LEADING EDGE: For explanations see Passaro et al. 2014 and 2018
        edgestart=1
        edgeend=1
        
        # --- Optional smoothing (introduced for SARAL with smooth = 3. Recommended values: s=3 for SARAL and Sentinel6_lrm, S=5 for ERS1/2)
        if self.smooth > 1 :
            wv0=D[index_originalbins]
            kernel_size = self.smooth
            kernel = np.ones(kernel_size) / kernel_size
            wv = np.convolve(wv0, kernel, mode='same')
        else:
            wv=D[index_originalbins] # old version of 'wv'
        
        Dwv=np.diff(wv)
        i = noisegates[0] #4 #Gate where the search starts (avoids wrap-up of ERS for example)
        
        while i<=np.size(wv)-5 :
            #In order to be the leading edge, it doesn't have to go below
            #the 1% of its maximum in the space of few gates
            if Dwv[i]>0.01 and wv[i+1]>thra and wv[i+2]> thra and wv[i+3]> thra and wv[i+4]> thra:
                edgestart=i # FOUND THE START OF THE LEADING EDGE
                
                while i<=np.size(Dwv)-4 :
                    #if the slope is negative for more than one gate, then the trailing edge
                    #is starting
                    if wv[i+1]>thrb:  # extra test added by FA for very large sea states and blooms (works OK with thrb>0.6 . 0.8 too large for some cases) 
                        edgeend=i # FOUND THE END OF THE LEADING EDGE    
                        if Dwv[i]<0 :
                            if Dwv[i+1]>0 and Dwv[i+2]>0 and Dwv[i+3]>0 : #if following gates are still growing, it's only a perturbation of the leading edge
                                i=i+1
                            else :
                                break
                        else :
                            i=i+1
                    else :
                        i=i+1
                break
            else :
                i=i+1
        
        leadextend=interpolation_factor
        gate2=index_originalbins[edgeend+leadextend]  
        #Gate2 is the end of the leading edge. One gate of tolerance is added for numerical reasons (for example if a leading-edge-only retracking is attempted)
        gate1=index_originalbins[edgestart]  #Gate1 is the start of the leading edge on the waveform

        self.gate2=gate2
        self.gate1=gate1

        start_a = time.time()
        

        if np.isnan(D[gate2])==0:

# STEP 5: RETRACKING WITH THE ALES RETRACKER
            
            x1_LESfive=np.empty(3)*np.nan  #Initialisation of vector of the estimated parameters (LESfive corresponds to the second retracking)
            x1_yang=np.empty(3)*np.nan     #Initialisation of vector of the estimated parameters (yang corresponds to the first retracking, restricted to the leading edge)
            Err_yang=1
    
            # STEP 5.1: FIRST PASS INITIAL CONDITIONS
            x_initial=(edgestart-1)*tau
            sigma_initial=(edgeend-edgestart)*tau/(2*np.sqrt(2))
            ampl_initial=2*np.mean(D[gate1:gate2])
            
            # STEP 5.2: FIRST PASS RETRACKING
            tol_Err_yang=Err_tolerance_vector
            x1_yang[2]=0 #initialising to start the while cycle
            growingdue=0 #variation of the number of gates considered
            SWH=0.1  # if not initialized here, SWH may not be defined with thrb=0.7  (maybe edgeend ins undefined ??) 
            while Err_yang>tol_Err_yang and (gate2+growingdue)<total_gate_number and sigma_initial<100 :
            #while cycle: check on the Fitting Error, check that the number of available gates is not being exceeding, check that initial condition has been defined as a number (sigma_inital<100, could be removed)    
                
                weightflag=0  #No weights are used in the first (leading-edge only) retracker
                this_weights=np.ones(np.shape(xdata))
                
                x1_yang, Wt_yang, exitflag_yang, Err, SWH =\
                        self.NM_fit( xdata[startgate:gate2+growingdue+1] , D[startgate:gate2+growingdue+1],\
                        self.xi*math.pi/180,self.tau,self.Theta,self.SigmaP,self.hsat,\
                        np.array([x_initial, sigma_initial, ampl_initial]),mission,this_weights[startgate:gate2+growingdue+1],weightflag,modelcost='brown_LS')
                
                #print('startgate',startgate)
                #print('gate2+growingdue+1',gate2+growingdue+1)
                #print('D[startgate:gate2+growingdue+1]',D[startgate:gate2+growingdue+1])
                #print('x1_yang',x1_yang)
                #print('SWH',SWH)
                #print('x_initial',x_initial)
                #print('exitflag_yang',exitflag_yang)



                if gate2>gate1+1 and gate1>0 and np.size(Wt_yang)>1 and gate1>startgate and gate2>startgate+1 :
                            # Fitting Error computed on the leading edge
                            Err_yang=np.sqrt( 1/float(gate2-gate1-1) * np.sum( (D[gate1+1:gate2]-Wt_yang[gate1-(startgate-1):gate2-startgate])**2 ))
                            
                self.Wt_all_yang[startgate:gate2+growingdue+1]=Wt_yang
                #If the fitted waveforms only contains zeroes,
                #it's not the good one...
                if np.sum(Wt_yang==0)==np.size(Wt_yang):
                    Err_yang=1
                
                
                if exitflag_yang==0: 
#                #if convergence is not reached, we might need more gates (growingdue+2) 
#                #We need a variable to define this in processing choices ... 
                    if mission.lower() == 'saral' or mission.lower() == 'saral_igdr' or mission.lower() == 'ers2_r_2cm' or mission.lower() == 'ers1':
                        if growingdue < 6 : #In the special case of Saral and ERS, where the search for a leading edge is based on a smoothed waveform, the convergence often fails 
                                            #for non-oceanic waveforms, therefore we limit the attempts to add more gates, in the interest of time
                            growingdue=growingdue+2
                        else:
                            break
                    else:
                        growingdue=growingdue+2
                    
                    if mission.lower() == 'sentinel6_lrm':
                        if growingdue < 6 : #In the special case of sentinel6, where the search for a leading edge is based on a smoothed waveform, the convergence often fails 
                                            #for non-oceanic waveforms, therefore we limit the attempts to add more gates, in the interest of time. Note that in WHALES for Ifremer
                                            #the same is also done for saral
                            growingdue=growingdue+2
                        else:
                            break
                    else:
                        growingdue=growingdue+2                        
                else: #else then we can stop the while cycle
                        break
                    

            self.D=D  #D is the normalised original waveform
  

                      
            del(growingdue)    
                    
            if np.isnan(x1_yang[2])==0 :
                Epoch_yang,Tau_yang,Sigma_yang,Au_yang=self.Conversion_NMbrown(x1_yang,tau,self.nominal_tracking_gate)
                SWH_yang=SWH #SWH after the first retracking 
                # self.SWH_yang=SWH_yang #SWH after the first retracking              
                Sigma0_yang=10*np.log10(Au_yang*normalize_factor) #De-normalise the amplitude and save it as db 
            else:
                    Epoch_yang=np.nan 
                    Tau_yang=np.nan
                    SWH_yang=np.nan
                    # self.SWH_yang=np.nan
                    Au_yang=np.nan
                    Sigma0_yang=np.nan 
                    Sigma_yang=np.nan
                             
                             
            end_a = time.time()
            #print(f"Time elapsed for section A: {end_a - start_a:.4f} seconds")                 
                             
                             
                             
             # STEP 5,3  WHALES EQUATION TO SET SECOND-PASS WINDOW: depending on the first estimation of SWH, a new width of the subwaveform is chosen, similarly 
                         #what is done in the ALES retracker (see Passaro et al. 2014 for the methodology). "stopgate" is the limit of the new subwaveform
             
            SWH_yang=np.real(SWH_yang)
            if SWH_yang<0:
                SWH_yang=0 # THE WHALES ALGORITHM IS DERIVED ONLY FOR POSITIVE SSB
            if SWH_yang< maxHs: #Estimations of high SWH in the ERS are too noisy, so we extend the retracking window for SWH_yang>10m instead of 15m
                gateafterLE=np.ceil(ALEScoeff0+ALEScoeff1*np.abs(SWH_yang))                  
                stopgate=np.ceil(x1_yang[0]/tau)+gateafterLE   
                #print('436')         
            else: #we consider maxHs m as the maximum possible SWH value: this used to be set to 15. 
                gateafterLE=np.ceil(ALEScoeff0+ALEScoeff1*maxHs)                   
                stopgate=np.ceil(x1_yang[0]/tau)+gateafterLE
                #print('440')
                #print(x1_yang[0])
            if interpolation_factor > 1 :
                if stopgate<=total_gate_number/interpolation_factor-1*interpolation_factor: #if we need more than the actual number of gates, we use the full waveform
                    pass
                else:
                    stopgate=total_gate_number/interpolation_factor-1*interpolation_factor #-1 because then we have the condition "(stopgate+growingdue)<total_gate_number"                
                    #print('444')
            else:
                if stopgate<=total_gate_number: #if we need more than the actual number of gates, we use the full waveform
                    pass
                else:
                    stopgate=total_gate_number-1 #-1 because then we have the condition "(stopgate+growingdue)<total_gate_number"  
                    #print('446')                  




            if xdata[-1]<stopgate*tau : #if the computed stopgate is bigger than the actual length of the vector
                                        #(can happen for very high SWH), then stopgate will be the end of the waveform
                stopgate=np.size(xdata)-1
            elif stopgate*tau<0: #if the computed stopgate is negative, something has failed in the fitting, we put the stopgate at the end of the waveform
                stopgate=np.size(xdata)-1
            else: # FIND ON THE OVERSAMPLED WAVEFORM THE POSITION CORRESPONDING TO THE COMPUTED STOPGATE                        
                #print(xdata)
                #print(stopgate*tau)
                stopgate=np.nonzero(xdata==stopgate*tau)[0]-1 #-1 because the first element of xdata is 0 and not 3.125



            #we don't want to stop before the end of the
            #leading edge
            
            if stopgate>=gate2:
                pass
            else:
                stopgate=gate2
                
            if np.isscalar(stopgate) == False: #Addition to check that the output of stopgate is a scalar
                stopgate=stopgate[0]

            
            
            
    
    
                   
             # STEP 5.4 SECOND PASS RETRACKING                  
            exitflag_LESfive=0
            Err_LESfive=1  
            tol_Err_LESfive=Err_tolerance_vector
                
            if Err_yang<1: #To make it faster: If the first pass was not able to retrack the leading edge, it's useless to go on
                    if Err_LESfive>tol_Err_LESfive:
                        growingdue=0 #variation of the number of gates considered
                        x1_LESfive[1]=0 #initializing to start the while cycle
                        while Err_LESfive>tol_Err_LESfive and (stopgate+growingdue)<total_gate_number and sigma_initial<100 : #first while cycle: check on the Fitting Error

                            weightflag=1
                            
                            #In WHALES, for each value of SWH there is a corresponding set of weights, which are defined only from the start to the end of the 
                            #leading edge of the synthetic waveforms generated in the Montecarlo simulation
                            

                            if self.weights_type == 0:
                                 this_weights=np.zeros(len(waveform))+1.0
                            elif self.weights_type == 1:
                                 weigths_SWH_vector=np.arange(0,10.5,0.5)
                                 select_weights=np.argmin(  np.abs(np.abs(SWH_yang)-weigths_SWH_vector)  )
                                #Select the right line of weights
                                 weights_select=self.weights[select_weights,:]
                                 index_nanweights=np.where(np.isnan(weights_select))[0]
                                 weights_select[index_nanweights]=self.weight_outsub #Transform the NaNs of the weight vector in ones
                                 index_startweight=np.where(self.weights_flag[select_weights,:]==1)[0] #identify the start and the end of the leading edge in the weight vector
                                 index_endweight=np.where(self.weights_flag[select_weights,:]==2)[0]
                                 index_startweight=index_startweight[0]*interpolation_factor #convert array to index
                                 index_endweight=index_endweight[0]*interpolation_factor #convert array to index

                                 if interpolation_factor > 1 :
                                     #Consider the possibility that the waveform is being oversampled and therefore the weights matrix has also to be oversampled
                                     # Get the original number of samples
                                     original_length = weights_select.shape[0]

                                     # Create a new index for the interpolated data
                                     original_indices = np.linspace(0, original_length - 1, num=original_length)
                                     new_indices = np.linspace(0, original_length - 1, num=total_gate_number)
       
                                     # Perform 1D interpolation
                                     weights_select = np.interp(new_indices, original_indices, weights_select)
                            
                                 this_weights=np.zeros(len(waveform))+self.weight_outsub
                                 if gate1+index_endweight-index_startweight<total_gate_number :
                                     this_weights[gate1:gate1+(index_endweight-index_startweight)]=weights_select[index_startweight:index_endweight]
                                     this_weights=1./this_weights

                                 #gatemax=np.min([total_gate_number-1,gate1+index_endweight-index_startweight])
                                 #this_weights[gate1:gatemax]=weights_select[index_startweight:index_startweight+gatemax-gate1]

                            elif self.weights_type == 2:
                                 x1_brown=x1_yang
                                 x1_brown[1]=max(x1_yang[1],minHs) 
                                 x1_brown[2]=1   # need to check the impact of this ... 
                                 this_weights=np.zeros(len(waveform))+self.weight_outsub
                                 Gamma=0.5 * (1/math.log(2))*np.sin(self.Theta)*np.sin(self.Theta) # antenna beamwidth parameter
                                 c=3.0*(10**8) #Light speed
                                 H=self.hsat
                                 Ri=6378.1363*(10**3) #Earth radius

                                 Gamma=0.5 * (1/math.log(2))*np.sin(self.Theta)*np.sin(self.Theta) # antenna beamwidth parameter
                                 Zeta=self.xi*math.pi/180
                                 b_xi = np.cos (2*Zeta) - ((np.sin(2*Zeta))**2)/Gamma
                                 a=( (4/Gamma)*(c/H) * 1/(1+H/Ri))
                                 c_xi=b_xi* ( (4/Gamma)*(c/H) * 1/(1+H/Ri))
                                 c_xi=c_xi/1000000000 #1/ns
                                 wbrown=0.001+wf_brown_eval(xdata,x1_brown,0.,Gamma,Zeta,c_xi,[1])/self.sqrtn # tau,PTR)
                                 this_weights[gate1:gate2]=wbrown[gate1:gate2]
                                 this_weights=1./this_weights


                            elif self.weights_type == 3:
                                 x1_brown=x1_yang
                                 x1_brown[2]=1   # need to check the impact of this ... 
                                 this_weights=np.zeros(len(waveform))+self.weight_outsub
                                 Gamma=0.5 * (1/math.log(2))*np.sin(self.Theta)*np.sin(self.Theta) # antenna beamwidth parameter
                                 c=3.0*(10**8) #Light speed
                                 H=self.hsat
                                 Ri=6378.1363*(10**3) #Earth radius

                                 Gamma=0.5 * (1/math.log(2))*np.sin(self.Theta)*np.sin(self.Theta) # antenna beamwidth parameter
                                 Zeta=self.xi*math.pi/180
                                 b_xi = np.cos (2*Zeta) - ((np.sin(2*Zeta))**2)/Gamma
                                 a=( (4/Gamma)*(c/H) * 1/(1+H/Ri))
                                 c_xi=b_xi* ( (4/Gamma)*(c/H) * 1/(1+H/Ri))
                                 c_xi=c_xi/1000000000 #1/ns
                                 wbrown=0.001+(wf_brown_eval(xdata,x1_brown,0.,Gamma,Zeta,c_xi,[1])/self.sqrtn)**2 # tau,PTR)
                                 this_weights[gate1:gate2]=wbrown[gate1:gate2]
                                 this_weights=1./this_weights
                            

                            # Launch the second retracking process
                            self.gate3=stopgate+growingdue+1
                            if self.costfunction == 'LS':
                                x1_LESfive, Wt_LESfive, exitflag_LESfive, Err, SWH =self.NM_fit( xdata[startgate:stopgate+growingdue+1] , \
                                D[startgate:stopgate+growingdue+1],self.xi*math.pi/180,tau,self.Theta,self.SigmaP,self.hsat,np.array([x_initial, sigma_initial, \
                                ampl_initial]),mission,this_weights[startgate:stopgate+growingdue+1],weightflag,modelcost='brown_LS')
                            else:
                                x1_LESfive, Wt_LESfive, exitflag_LESfive, Err, SWH =self.NM_fit( xdata[startgate:stopgate+growingdue+1] , \
                                D[startgate:stopgate+growingdue+1],self.xi*math.pi/180,self.tau,self.Theta,self.SigmaP,self.hsat,np.array([x_initial, sigma_initial, \
                                ampl_initial]),mission,this_weights[startgate:stopgate+growingdue+1],weightflag,modelcost='brown_ML')
                            
                                     
                            self.Wt_all_LESfive[startgate:stopgate+growingdue+1]=Wt_LESfive  #This is the fitted subwaveform
        
                                
                            if gate2>gate1+1 and gate1>0 and np.size(Wt_LESfive)>1 and gate1>startgate and gate2>startgate+1 :
                                 #Fitting error on the leading edge
                                Err_LESfive=np.sqrt( 1/float(gate2-gate1-1) * np.sum( (D[gate1+1:gate2]-Wt_LESfive[gate1-(startgate-1):gate2-startgate])**2 ))
                                

                                
                                
                            if exitflag_LESfive==0 and tol_Err_LESfive==Err_tolerance_vector : #if convergence is not reached, we might need more gates (growingdue+1), otherwise we prepare to terminate the while cycle
        
                                growingdue=growingdue+2

                                    
                                    
                            else: # then we can stop the while cycle
                                growingdue=total_gate_number-gate2;
                                
                                
                                
                                 
                            #END OF THE SECOND WHILE CYCLE

            end_b = time.time()
            #print(f"Time elapsed for section B: {end_b - end_a:.4f} seconds") 

                 
                
            if np.isnan(x1_LESfive[2])==0 :
                Epoch_LESfive,Tau_LESfive,Sigma_LESfive,Au_LESfive=self.Conversion_NMbrown(x1_LESfive,tau,self.nominal_tracking_gate)
                SWH_LESfive=SWH
    
                Sigma0_LESfive=10*np.log10(Au_LESfive*normalize_factor)  #db
            else :
                Epoch_LESfive=np.nan 
                Tau_LESfive=np.nan
                SWH_LESfive=np.nan
                Au_LESfive=np.nan
                Sigma0_LESfive=np.nan 
                Sigma_LESfive=np.nan
            
            
            

#
# STEP 5.5 FINAL OUTPUT
#
        if np.isnan(Epoch_LESfive)==0:
            
            self.Epoch=Epoch_LESfive
            self.SWH=SWH_LESfive  #in m
            self.SWH_yang=SWH_yang
            self.Sigma=Sigma_LESfive #Rising time of the leading edge in ns
            self.Amplitude=Sigma0_LESfive#  Instead of the normalised amplitude(Au_LESfive), we are outputting the backscatter coefficient (before corrections)
            #self.Norm_Amplitude=Au_LESfive
            #self.leading_edge = np.nan
            self.Error = Err_LESfive  #Fitting Error on the leading edge
            
            self.range = self.uralt  + self.Epoch   #Range in m  
            self.D=D
            self.this_weights=this_weights
        else:
            self.range = np.nan
            self.uncssh = np.nan            
            self.SWH=np.nan 
            self.SWH_yang=np.nan
            self.Sigma=np.nan 
            self.Amplitude=np.nan 
            self.Norm_Amplitude=np.nan 
            #self.leading_edge = np.nan
            self.Error = np.nan 
            self.D=D
            self.this_weights=D*0+np.nan
            self.gate3=total_gate_number

        self.model = self.Wt_all_LESfive.copy()  #*normalize_factor  #Fitted Waveform
            
        
       

        
        
        
        
        
        
        
        
        
        

        
        
