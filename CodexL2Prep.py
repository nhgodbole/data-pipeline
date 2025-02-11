#!/usr/bin/env python

"""
CODEX L2 Prep - Test version.
Created on Mon Jan  6 19:14:14 2025

@author:
    mcasti
    The Catholic University of America at NASA/GSFC
    email: marta.casti@nasa.gov
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from astropy.io import fits

from codex_pipeline import *


def get_dns(data):
    # Get exposure time.  #TO DO: use exposure time in seconds
    exp_time = data['hdr']['DIT']/1000
    # Divide image for the exposure time
    data_dns = {}
    data_dns['hdr'] = data['hdr']
    data_dns['img'] = data['img']/exp_time
    return(data_dns)


def apply_vign(vign, data):
    # APPLY VIGNETTING FUNCTION
    data_novign = {}
    pols = ['pol_1', 'pol_2', 'pol_3', 'pol_4']
    for pol in pols:
        vign[pol][np.where(vign[pol] == 0)] = 1 * 10 ** -10
        data_novign[pol] = data[pol]*(1/vign[pol])
    return(data_novign)


if __name__ == '__main__':
    # Get L1 data product to be processed from its directory
    input_dir = 'Input_L2'
    output_dir = 'Output_L2'
    ########################
    # Load CODEX parameters
    ########################

    # Darks
    dark_dir = get_path('Darks')

    # Detector characteristic
    detector_feat = get_param('detector_feat')
    # Load Demodulation Tensor
    demod_path = get_path('dem_tens')
    demod = get_demod(demod_path)

    # Load Vignetting
    vign_path = get_path('vign')
    vign = get_vign(vign_path)
    # Upload QE data as Pandas dataframe
    path_qe = get_path('qe')
    wavelength = 423.5
    qe_filter = get_qe(path_qe, wavelength)  # TO DO: take filter info from hdr

    # Get CODEX image to be processed form Input folder
    data = get_img_from_dir(input_dir)
    filename = list(data.keys())[0]
    corona_img = data[filename]
    texp = corona_img['hdr']['DIT']

    ######################
    #    DARK REMOVAL    #
    ######################
    # Get available images of dark for the same
    dark_img = get_dark(texp, dark_dir)
    # Remove dark from the input image
    corona_nodark = remove_dark(corona_img, dark_img)
    # Flip image up-down
    corona_nodark['img'] = np.flipud(corona_nodark['img'])

    ####################
    #    T_EXP NORM    #
    ####################
    # Divide for the exposure time
    corona_nodark_dns = get_dns(corona_nodark)
    # Separate the 4 polarized images. NOTE: the four images should be already
    # separated within the fits file - TO DO: change accordingly
    corona_nodark_dns = get_pol_imgs(corona_nodark_dns)

    ####################
    #    VIGNETTING    #
    ####################
    # Apply vignetting
    corona_novign = apply_vign(vign, corona_nodark_dns)

    ######################
    #    POL ANALYSIS    #
    ######################
    # Polarization analysis on the images
    plot_flag = 0
    stokes, pol_aolp, pol_dolp = analyze_pol(corona_novign, demod,
                                             parameters['detector_pol'],
                                             parameters['mask_pol'],
                                             detector_feat['gain'],
                                             qe_filter, plot_flag, output_dir)
    # Derive polarized brightness image from Stokes vector elements
    pB_img = get_pb(stokes)
    # Plot polarized brightness image
    plot_pB(pB_img)
