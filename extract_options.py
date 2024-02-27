# -*- coding: utf-8 -*-
"""
Created on November 2018

The code launches the WHALES retracker using original mission files

Works for the following missions: 

@author: Marcello Passaro
"""
import sys
sys.path.append('/home/mdecarlo/Documents/PROJETS/codes_Python/whales/whales/')
import argparse
# import cmath
import netCDF4
from netCDF4 import Dataset
# import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
# import matplotlib.pyplot as plt
import scipy.io
import os
# import matplotlib.animation as manimation
import time
from compute_instr_corr_SWH_WHALES import compute_instr_corr_SWH_WHALES
# import sys
# from read_functions import wf_reader

from Retracker_MP import *

from WHALES_withRangeAndEpoch import *

from scipy.io import matlab
# import pandas as pd

from import_weights_mat import import_weights_mat


def get_options():
    parser = argparse.ArgumentParser(
        description='Retrack a SGDR file with WHALES')

    parser.add_argument(
        '-m', '--mission', type=str,
        choices=['envisat', 'jason1', 'jason2', 'jason3', 'saral', 'cs2_lrm',
                 'jason3f'],
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
    return parser.parse_args()


options = get_options()

np.savez('options',options)
