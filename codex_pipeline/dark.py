"""
CODEX Dark images functions.
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

def get_dark(exp_time_img, dark_dir):

    dark_imgs_metadata = pd.read_csv(os.path.join(dark_dir, 'metadata.csv'))
    return _get_dark(exp_time_img, dark_imgs_metadata, dark_dir)

def _get_dark(exp_time_img, dark_metadata, dark_dir):
    '''PURPOSE: Look or generate dark image with the same exposure time of the
            provided image.
        INPUT:
            exp_time_img (float) = query criteria: exposure time of the image
                to subtract dark signal.
            dark_metadata (dict) = contais the metadata of the dark images.
            dark_dir (str) = directory where dark images are stored.
        OUTPUT:
            dark_filenames (list) = filename of the dark image to be used to
                subtract dark.
        INVOKED FUNCTIONS:
            - util.search_imgs
            - UTIL.get_img_from_list
            - get_avgd_img
            - interp_imgs'''
    # Get metadata of dark images acquired with the same exposure time of the
    # image
    same_exp_time_dark = search_imgs(dark_metadata, 'DIT', 'equal',
                                     exp_time_img)
    # Get the number of available dark images
    n_darks = len(same_exp_time_dark.index)
    if n_darks == 1:
        # Only one dark image with the needed exposure time is available -->
        # return the image
        dark_filenames = list(same_exp_time_dark.loc[:, 'FILENAME'])
        # Collect the matching dark image from the repository
        dark_img = get_img_from_list(dark_dir, dark_filenames)
    elif n_darks > 1:
        # More than 1 dark image with the needed are available
        dark_filenames = list(same_exp_time_dark.loc[:, 'FILENAME'])
        # Collect the matching dark images from their repository
        dark_imgs = get_img_from_list(dark_dir, dark_filenames)
        # Average all dark frame with the same exposure time to generate a
        # unique dark image
        dark_img = get_avgd_img(dark_imgs, same_exp_time_dark, 'dark')
    else:
        # No dark images with the needed exposure time are available
        exp_time_1, img_1, exp_time_2, img_2 = get_nearest_imgs(dark_metadata,
                                                                exp_time_img,
                                                                dark_dir)
        dark_img = {}
        dark_img['img'] = interp_imgs(img_1['img'], exp_time_1, img_2['img'],
                                      exp_time_2, exp_time_img)
        dark_img['hdr'] = fits.Header()
        # raise ValueError('NO dark image found for this exposure time')
    return(dark_img)


def get_nearest_imgs(dark_metadata, exp_time_img, dark_dir):
    ''' PURPOSE: Select the images nearest in time to the analyzed one
        INPUT:
            - dark_metadata (dataframe) = contains metadata of the available
                dark images
            - exp_time_img (float) = exposure time to aim at
            - dark_dir (str) = directory containing the images
        OUTPUT:
            - exp_time_1 (float) = exposure time of the first image to
                interpolate
            - dark_img_1 (2D array) = first image nearest in time to the
                provided one
            - exp_time_2 (float) = exposure time of the second image to
                interpolate
            - dark_img_2 (2D array) = second image nearest in time to the
                provided one
        INVOKED FUNCTIONS:
            - find_nearest
            - pick_dark_img'''
    # Collect exp. time of the available images
    exp_times = np.unique(list(dark_metadata.loc[:, 'DIT']))
    # Find exposure time of the dark images with the exposure time closer to
    # the image to analyze
    found_exp_times = find_nearest(exp_times, exp_time_img)
    # Sort found exposure times
    found_exp_times.sort()
    # First exposure time: find all images acquired with that specific exposure
    # time
    exp_time_1 = found_exp_times[0]
    dark_img_1 = pick_dark_img(exp_time_1, dark_metadata, dark_dir)
    # Second exposure time: find all images acquired with that specific exp
    # time
    exp_time_2 = found_exp_times[1]
    dark_img_2 = pick_dark_img(exp_time_2, dark_metadata, dark_dir)
    return(exp_time_1, dark_img_1, exp_time_2, dark_img_2)


def find_nearest(array, value):
    ''' PURPOSE: Find the two elements in array closest to the assigned value
        INPUT:
            - array (2D array) = contains the element among wich select the two
                nearest to value
            - value (float) = assign reference value
        OUTPUT:
            - value1 (float) = first element closest in value to the reference
                one
            - value2 (float) = second element closest in value to the reference
                one
        INVOKED FUNCTIONS:
            None'''
    # TO DO: add checks for possible failures
    array = np.asarray(array)
    # Look for the index of the first-nearest element
    idx1 = (np.abs(array - value)).argmin()
    value1 = array[idx1]
    # Remove the first-nearest element from the array
    array = np.delete(array, idx1)
    # Look for the index of the second-nearest element
    idx2 = (np.abs(array - value)).argmin()
    value2 = array[idx2]
    return([value1, value2])


def get_avgd_img(data, metadata, feat):
    ''' PURPOSE: Compute the averaged image of the images provided as input.
        INPUT:
            - data (dict) = contains images and the related metadata
            - metadata (dataframe) = metadata related to the images in data.
                These are the UPDATED metadata
            - feat (string) = to be reported in the filename of the avg image
        OUTPUT:
            - avg_img (dict) = contains the image resulting from the average
                and its metadata. Structure is
                    avg_img{ [filename] {['img'], ['hdr']} }
        INVOKED FUNCTIONS:
            - avg_imgs '''
    # Initialized output dictionary
    avg_img = {}
    # IMAGE
    # Evaluate the average of the images and the associated standard deviation
    avg_img['img'], avg_img['std'] = avg_imgs(data)
    # HEADER
    # Generate new header for the resulting image starting from the metadata of
    # the first considered image
    new_header = metadata.iloc[0].to_dict()
    # Remove first key of the dictionary, since it is the pandas dataframe
    # index (TO DO: Verify)
    first_key = list(new_header.keys())[0]
    del(new_header[first_key])
    # Generate a filename for the new image
    new_header['FILENAME'] = feat+'_avg'
    # Create 'HISTORY' keyword
    new_header['HISTORY'] = 'Result of the average of images'
    # TO DO: add the acquisition time of the first and the last images used to
    # create the averaged image for n in np.arange(1,6,1):
    #    new_header['TEMP'+str(n)]=np.mean(metadata.loc[:,'TEMP'+str(n)])
    # Save new header in the dictionary
    avg_img['hdr'] = new_header
    return(avg_img)


def pick_dark_img(exp_time, dark_metadata, dark_dir):
    ''' PURPOSE:
        INPUT:
        OUTPUT:
        INVOKED FUNCTIONS:
            - search_imgs
            - util.get_img_from_list
            - util.get_avgd_img'''
    exp_time_imgs = search_imgs(dark_metadata, 'DIT', 'equal', exp_time)
    n_darks_1 = len(exp_time_imgs.index)
    if n_darks_1 == 1:
        # Only one dark image with the needed exposure time is available -->
        # return the image
        dark_filenames = list(exp_time_imgs.loc[:, 'FILENAME'])
        # Collect the matching dark image from the repository
        dark_img = get_img_from_list(dark_dir, dark_filenames)
    else:  # There are images, more than one
        # More than 1 dark image with the needed are available
        dark_filenames = list(exp_time_imgs.loc[:, 'FILENAME'])
        # Collect the matching dark images from their repository
        dark_imgs = get_img_from_list(dark_dir, dark_filenames)
        # Average all dark frame with the same exposure time to generate a
        # unique dark image
        dark_img = get_avgd_img(dark_imgs, exp_time_imgs, 'dark')
    return(dark_img)


def interp_imgs(img_1, exp_time_1, img_2, exp_time_2, exp_time):
    ''' PURPOSE: Interpolate images to create a syntetic image with an exposure
            time in between
        INPUT:
            - img_1 (2D array) = first image nearest in time to the provided
                one
            - exp_time_1 (float) = exposure time of the first image to
                interpolate
            - dark_img_2 (2D array) = second image nearest in time to the
                provided one
            - exp_time_2 (float) = exposure time of the second image to
                interpolate
        OUTPUT:
            - img (2D array) = syntetic image resulting from the interpolation
                for the considered exposure time
        INVOKED FUNCTION: '''
    aDiff1D = np.reshape(img_2 - img_1, (1, np.size(img_1)))
    aDiffRepeated = np.repeat(aDiff1D, np.size(exp_time), axis=0)
    aStep = aDiffRepeated*((exp_time-exp_time_1)/(exp_time_2-exp_time_1))
    aInterp = np.repeat(np.reshape(img_1, (1, np.size(img_1))),
                        np.size(exp_time), axis=0)+aStep
    aInterp = np.reshape(aInterp, (np.size(exp_time), img_1.shape[0],
                                   img_1.shape[1]))
    img = aInterp[0]
    return(img)


def avg_imgs(data):
    ''' PURPOSE: Generate an image which is the average of the images provided
            in input.
        INPUT:
            - data (dict) = contains images and the related metadata.
                Structure is the following:
                    data = { [filename1.fits] = { ['img']=[…], ['hdr']=[…] },
                        [filename2.fits] = {…} }
        OUTPUT:
            avgd_img (2D array) = Average of the images provided as input.
            std_img (2D array) = Standard deviation of the averaged image.
        INVOKED FUNCTIONS:
            None'''
    imgs = {}
    # Extract images from dictionary
    for filename in list(data.keys()):
        imgs[filename] = data[filename]['img']
    # Stack images contained in the data dictionary to compute the mean
    stacked_imgs = np.stack(list(imgs.values()))
    # Average the received images and save the result in the output dictionary
    avgd_img = stats.trim_mean(stacked_imgs, proportiontocut=0.1, axis=0)
    # Standard deviation for the averaged image
    std_img = np.std(stacked_imgs, axis=0)
    return(avgd_img, std_img)


def remove_dark(data, data_dark):
    '''PURPOSE: Remove background and dark from images.
        INPUT:
            - data (dict) = contains the image from which remove the dark
                signal. The structure is the following:
                    data {['img']:[], ['hdr']:[]}
            - data_dark (dict) = contains the dark image and its header.
                The structure is the following:
                    data_dark {['img']:[], ['hdr']:[]}
        OUTPUT:
            - no_dark (dict) = contains the image without dark and its updated
                metadata.The structure is the following:
                    no_dark {['img']:[], ['hdr']:[]}
        INVOKED FUNCTIONS:
            - subtract_imgs'''
    img = data['img']
    img_hdr = data['hdr']
    dark = data_dark['img']
    # Initialize output vector
    no_dark = {}
    # IMAGE
    # Change type to array, subtract dark and equal to zero possible negative
    # pixels
    nodark_img = subtract_imgs(img, dark)
    # HEADER
    # Copy header from the image in input and modify accordingly
    nodark_hdr = fits.Header()
    nodark_hdr = img_hdr
    # Add History line describing the process to the data product metadata
    nodark_hdr['HISTORY'] = 'Result of dark subtraction'
    # nodark_hdr['FLATS'] = img_hdr['FILENAME']
    # Add list of images used to generate the data image
    # nodark_hdr['DARKS'] = dark_hdr['FILENAME']
    # Add list of dark images used to generate the dark image
    # nodark_hdr['FILENAME'] = 'codex_'+level+'_'+datetime.utcnow().strftime('%Y%m%dT%H%M%S')+'.fits'   #TO DO: generalize 'cal'
    # nodark_hdr['DATE'] = datetime.utcnow().strftime('%Y%m%dT%H%M%S')
#     if nodark_hdr['HISTORY'] :
#         nodark_hdr['HISTORY1'] = 'Result of dark subtraction'
#     else:
#         nodark_hdr['HISTORY'] = 'Result of dark subtraction'
    # Fill ouput dictionary
    no_dark['img'] = nodark_img
    no_dark['hdr'] = nodark_hdr
    return(no_dark)


def subtract_imgs(img1, img2):
    '''PURPOSE: Subtract img2 from img1.
        INPUT:
            img1 (2D array) = first image
            img2 (2D array) = second image to be subtracted from img1
        OUTPUT:
            new_img (2D array) = result of the image subtraction (img1-img2)
        INVOKED FUNCTIONS:
            None'''
    new_img = np.subtract(np.array(img1, dtype=np.float32),
                          np.array(img2, dtype=np.float32))
    new_img[new_img < 0] = 0
    return(new_img)