# seastatecci_whales

The WHALES retracker is an algorithm aimed at estimating the significant wave height from waveforms acquired by satellite altimeters.

The code is an implementation of the algorithm described in:

Passaro M. et al., 2021: Algorithm Theoretical Basis Document (ATBD), Sea State Climate Change Initiative, European Space Agency, accessible from https://climate.esa.int/media/documents/Sea_State_cci_ATBD_v3.0-signed.pdf"

## Credits: 
- Marcello Passaro (DGFI-TUM) is the author of the WHALES algorithm, with support of Paolo Cipollini (previously at National Oceanography Centre, UK, now at the European Space Agency) and Fabrice Ardhuin (LOPS/Ifremer, Brest). The WHALES algorithm is an evolution of ALES (Passaro et al., 2014)

- Passaro M., Cipollini P., Vignudelli S., Quartly G., Snaith H.: ALES: A multi-mission subwaveform retracker for coastal and open ocean altimetry. Remote Sensing of Environment 145, 173-189, 10.1016/j.rse.2014.02.008, 2014

- Passaro, M., Dodet, G., Ardhuin, F., & Cipollini, P. (2025). WHALES: an optimized retracker for satellite radar altimeter waveforms in sea state applications. arXiv preprint arXiv:2505.12881.



## Getting started

First, clone this repository:

    $ git clone https://github.com/ne62rut/whales

Go into newly created directory, to which the repository was cloned

    $ cd seastatecci_whales

Then, create a virtual environment using conda and taking as reference the environment file provided:

    $ environment.yml

(This should install all the dependencies/packages/python-version (including version numbers) that are set in environment.yml)

Switch to into the "isolated" virtual environment 

    $ conda activate seastatecci_whales


## Usage and command line options

The retracker shall be launched from the command line through the file python_WHALES_launcher.py with the following suggested options (sample files):

$ python python_WHALES_launcher.py -m MISSION -s 5 -i E2_REAP_ERS_ALT_2S_20000718T222234_20000718T235223_RP01.NC -o /output/

Supported missions are: ers1, ers2, sentinel6_lrm, swot, saral, cfosat, jason1, jason2, jason3f 

The parameter -s defines the size of the spatial smoothing of the waveform needed to correctly detect the leading edge. The suggested options are -s 5 for ers1 and ers2; -s 3 for sentinel6_lrm, swot, saral, cfosat. For the other missions no smoothing has been tested as yet, therefore the suggestion is not to add this parameter.

The retracker code is contained in:

    $ WHALES_withRangeAndEpoch.py
    
For Envisat, Saral, and Cryosat-2 NO INSTRUMENTAL CORRECTION IS USED.
For Jason-3, an instrumental correction to be ADDED to the retracked significant wave height can be added using the function 

    $ compute_instr_corr_SWH_WHALES.py

and an external correction model that associate each value of SWH to the correction, as in:

    $ SWHinstrcorr_WHALES_jason3SGDRd.mat
    
Note that the instrumental correction is simply based on the comparison with the one applied in the MLE3 retracker of the standard product. 
    
    
On the display, you will see a waveform counter for each successful retrack.
The launcher will save a NetCDF file with the same name of the original product.

### Optional smoothing (recommended for ERS1/2, SARAL ... and high wave height values) 
-s : forces smoothing before leading edge detection (see Passaro et al. 2025 for details)

-S : smoothing before leading edge detection but only for large wave heights
The definition of the leading edge is sensitive to noise and can produce leading edges that are too short and stop before the proper end of the leading edge. Smoothing the waveform reduces that 
sensitivity. Note that the -S option does two things (that are specific for large wave heights): 
- smoothing the waveform
- using theoretical weights (otherwise obtained with -w 2): this adjusts the weights beyond SWH=10 which would have been used otherwise (see Passaro et al. 2025) 

Here are some examples: with smoothing for Jason-3 but only for SWH > 8 m: 

    $ python my_folder/src/python_WHALES_launcher.py -m jason3f -i JA3_GPS_2PfP180_122_20210101_012654_20210101_022307.nc  -S 3 -o TEST_v5b


### Output for debug: leading edge indices and weights
If you want to get more information (for debug or other purposes) in the NetCDF output: the -d 1 option will add leading edge indices and weights and the waveform in the NetCDF output file
    

## Disclaimer

The only parameter provided for which the author is responsible is the Significant Wave Height.

Range and Backscatter Coefficients are provided without any verification.

The Backscatter Coefficient from WHALES is noisier than the Backscatter Coefficient from a full-waveform retracker, since WHALES subwaveform strategy is NOT adapted to the estimation of the backscatter.

The backscatter coefficient is the sum of the retracked signal amplitude, a scaling factor and an atmospheric correction. IN CRYOSAT-2, the backscatter coefficient provided corresponds ONLY to the retracked signal amplitude (of the waveform pwr_waveform_20_ku from the L1B files), since the scaling factor and the atmospheric correction are not provided in the L1B product used to extract the waveforms.

The WHALES Range of Cryosat-2 is provided using as on-board tracker the window delay field (window_del_20_ku). To the best of the author's knowledge, instrumental corrections that are usually applied to the on-board tracker of the other missions are NOT applied in the window delay field. These instrumental corrections are not available in the L1B files of Cryosat-2. 

## Known Issues







