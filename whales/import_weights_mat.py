# -*- coding: utf-8 -*-
"""


@author: Marcello Passaro
"""


import h5py  
import numpy as np


def import_weights_mat(my_path_weights)     :

    mat_weights = h5py.File(my_path_weights,'r')
    residual_std=np.transpose(mat_weights['residual_tot'].value) 
    flag_edges=np.transpose(mat_weights['flag_edges'].value   )
    
    return residual_std, flag_edges
