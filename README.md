# seastatecci_whales

A project for the editing and testing of the WHALES retracker

## Getting started

First, clone this repository:

    $ git clone https://gitlab.lrz.de/ne62rut/seastatecci_whales.git

Go into newly created directory, to which the repository was cloned

    $ cd seastatecci_whales

Then, create a virtual environment using conda and taking as reference the environment file provided:

    $ environment.yml

(This should install all the dependencies/packages/python-version (including version numbers) that are set in environment.yml)

Switch to into the "isolated" virtual environment 

    $ conda activate seastatecci_whales


## Usage

To run the retracker, open 

    $ python_WHALES_launcher.py

Before launching it, select the following parameters: saving_directory (default is the same directory where the launcher is contained), saving_name (please do not add file extension), filename (original file to be retracked) and mission (choose between envisat, jason1, jason2, jason3, saral, cs2_lrm).

The retracker code is contained in:

    $ WHALES_withRangeAndEpoch.py
    
For Envisat, Saral and Cryosat-2 NO INSTRUMENTAL CORRECTION IS USED.
For Jason-3, an instrumental correction to be ADDED to the retracked significant wave height can be added using the function 

    $ compute_instr_corr_SWH_WHALES.py

and an external correction model that associate each value of SWH to the correction, as in:

    $ SWHinstrcorr_WHALES_jason3SGDRd.mat
    
Note that the instrumental correction is simply based on the comparison with the one applied in the MLE3 retracker of the standard product. A new correction had to be computed by PML, but the performances were worse (according to the Round Robin results). 
    
    
    
On the display, you will see a waveform counter for each successful retrack ("Retracking waveform XXXX of YYYY"). Please note that the numbers displayed on the counter are not correct.
The launcher will save a NetCDF file of the kind produced for the Round Robin. It can be tested to work with the following test file from Jason-3:

    $ JA3_GPS_2PdP054_149_20170801_175030_20170801_184643.nc
    

## Disclaimer

The only parameter provided for which the author is responsible is the Significant Wave Height. Nevertheless, the only validation performed by the author is relative to Jason-3 in the framework of the Round Robin of the Sea State CCI. For the other missions, only visual test checks have been performed. The author is fully available to correct possible problems based on the feedback of the Validation Team of the Sea State CCI.

Range and Backscatter Coefficients are provided without any verification.

The Backscatter Coefficient from WHALES is noisier than the Backscatter Coefficient from a full-waveform retracker, since WHALES subwaveform strategy is NOT adapted to the estimation of the backscatter.

The backscatter coefficient is the sum of the retracked signal amplitude, a scaling factor and an atmospheric correction. IN CRYOSAT-2, the backscatter coefficient provided corresponds ONLY to the retracked signal amplitude (of the waveform pwr_waveform_20_ku from the L1B files), since the scaling factor and the atmospheric correction are not provided in the L1B product used to extract the waveforms.

The WHALES Range of Cryosat-2 is provided using as on-board tracker the window delay field (window_del_20_ku). To the best of the author's knowledge, instrumental corrections that are usually applied to the on-board tracker of the other missions are NOT applied in the window delay field. These instrumental corrections are not available in the L1B files of Cryosat-2. 

## Known Issues







