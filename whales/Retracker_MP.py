import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import traceback

from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.integrate import quad

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

    ' Additional Retracker Output from ALES '
    gate1 = None
    gate2 = None
    Wt_all_yang = None
    Wt_all_LESfive = None
    D = None
    Epoch = None
    SWH = None
    Amplitude = None
    Norm_Amplitude = None
    uncssh = None
    classifier = None
    interpolator = None
    Sigma = None

    ' Speed of light [m/s] '
    c = 299792458

    def __init__(self, config):
        if 'waveform' in config:
            self.waveform = np.array(config['waveform']).astype(float).tolist()

        if 'waveform_last' in config:
            self.waveform_last = np.array(config['waveform_last']).astype(float).tolist()

        if 'waveform_current' in config:
            self.waveform_current = np.array(config['waveform_current']).astype(float).tolist()

        if 'waveform_next' in config:
            self.waveform_next = np.array(config['waveform_next']).astype(float).tolist()

        if 'refbin' in config:
            self.refbin = config['refbin']
        if 'binsize' in config:
            self.binsize = config['binsize']
        if 'uralt' in config:
            self.uralt = config['uralt']
        if 'threshold' in config:
            if 0.0 <= float(config['threshold']) <= 100.0:
                self.threshold = float(config['threshold']) / 100.0
            else:
                print("Error: Threshold must be between 0 and 100%")
                sys.exit(0)
        if 'scale' in config and config['scale'] == 1:
            self.scale = 1

        ' Improved threshold '
        if 'ref_range' in config:
            self.ref_range = config['ref_range']
        if 'option' in config:
            self.option = config['option']

        if 'factor' in config:
            self.factor = config['factor']

        # MP added
        ' ales_rev'
        if 'mission' in config:
            self.mission = config['mission']
        if 'hsat' in config:
            self.hsat = config['hsat']
        if 'theta0' in config:
            self.theta0 = config['theta0']
        if 'xi' in config:
            self.xi = config['xi']
        if 'doppler' in config:
            self.doppler = config['doppler']
        if 'Epoch' in config:
            self.Epoch = config['Epoch']
        if 'SWH' in config:
            self.SWH = config['SWH']
        if 'Amplitude_initial' in config:
            self.Amplitude_initial = config['Amplitude_initial']
        if 'classifier' in config:
            self.classifier = config['classifier']
        if 'weights' in config:
            self.weights = config['weights']
        if 'weights_flag' in config:
            self.weights_flag = config['weights_flag']

        ' debug flag '
        if 'debug' in config:
            self.debug = int(config['debug'])

        if self.__class__.__name__ == "Threshold" and self.threshold is None:
            print("Error: Threshold retracker needs argument 'threshold' [0, ..., 100]")
            sys.exit(0)
        if self.__class__.__name__ == "ImprovedThreshold" and self.threshold is None:
            print("Error: ImprovedThreshold retracker needs argument 'threshold' [0, ..., 100]")
            sys.exit(0)

    def plot(self):
        x_data = list(range(len(self.waveform)))
        y_data = self.waveform
        plt.plot(x_data, y_data, 'b.-')