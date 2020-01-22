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

The retracker code is contained in:

    $ WHALES_withRangeAndEpoch.py

Additionally, an instrumental correction to be ADDED to the retracked significant wave height can be added using the function 

    $ compute_instr_corr_SWH_WHALES.py

and an external correction model that associate each value of SWH to the correction, as in:

    $ SWHinstrcorr_WHALES_jason3SGDRd.mat
    
Note that the instrumental correction is simply based on the comparison with the one applied in the MLE3 retracker of the standard product. A new correction had to be computed by PML, but the performances were worse (according to the Round Robin results). Also note that Envisat does not have an instrumental correction. 
    
To run the retracker, launch 

    $ python_WHALES_launcher.py
    
    
On the display, you will see a waveform counter for each successful retrack ("Retracking waveform XXXX of YYYY"). Please note that the numbers displayed on the counter are not correct.
The launcher will save a NetCDF file of the kind produced for the Round Robin. It is already set to work with the following test file from Jason-3:

    $ JA3_GPS_2PdP054_149_20170801_175030_20170801_184643.nc





