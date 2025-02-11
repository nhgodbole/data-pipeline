"""
CODEX utility functions.
Created on Mon Jan  6 19:14:14 2025

@author:
    mcasti
    The Catholic University of America at NASA/GSFC
    email: marta.casti@nasa.gov
"""


import os, csv
import functools
from collections import Counter
from os.path import dirname, isfile, join as pjoin
from os import listdir
from astropy.io import fits
from astropy.visualization import hist
from astropy.stats import mad_std
from matplotlib import pyplot as plt
from matplotlib import patches as patches
from datetime import datetime
from scipy import stats
import numpy as np
import pandas as pd

# __all__ = ["get_img_from_list", "get_img_from_dir", "search_imgs", "query_imgs", "avg_imgs", 
#            "imgs_stats", "subtract_imgs", "gen_datacube", "collect_metadata", "save_fits",
#           "get_pol_imgs", "gen_stats", "get_stats"] 


def get_img_from_list(file_dir, filenames):
    ''' PURPOSE: Collect in a dictionary the  file listed in filenames and stored in file_dir (directory). 
            The dictionary contains the image and its header.
        INPUT: 
            file_dir (str) = path pointing to the directory where the fits file to be read are saved.
            filenames (list) = list of files to be read in the directory.
        OUTPUT:
            files (dict) = { [filename1.fits] = { ['hdr']=[…], ['img']=[…] }, 
                        [filename2.fits] = {…}; ... }
            flag (bin) = Flag indicating if the process was successfully executed
        INVOKED FUNCTIONS: 
            None. '''
    # Check the list is not empty
    if len(filenames)>0:
        # Initialize output dictionary
        files = {}
        # Loop on filenames contained in the list
        for filename in filenames:
            files[filename] = {}
            # Open fits file, verify and fix header
            hdul = fits.open(os.path.join(file_dir, filename), checksum=True)
            hdul[0].verify('fix')
            files[filename]['hdr'] = hdul[0].header
            files[filename]['img'] = hdul[0].data 
            # Assign value to validity flag
            # flag = 1
    else:
        # Assign value to validity flag
        # flag = 0
        # Rise error if the directory in input does not contain fits files 
        raise ValueError('Empty file list')    
    return(files)


def get_img_from_dir(file_dir):
    ''' PURPOSE: Create a list of all fits files contained within the assigned directory and 
            generate a dictionary containing, for each fits file in the directory, the header
            and the 2D array rappresenting the image.
        INPUT: 
            file_dir (str) = path pointing to the directory where the fits file are stored.
        OUTPUT: 
            files (dict) = nested dictionary, contains all read fits file as header and 2D array. The structure is the following: 
                files = { [filename1.fits] = { ['hdr']=[…], ['img']=[…] }, 
                            [filename2.fits] = {…}; ... }
        INVOKED FUNCTIONS:
            - get_img_from_list '''
    # Look for fits file in the assigned directory and collect their filename as a list
    filenames = [file for file in listdir(file_dir) if file.endswith(".fits")]
    # Create a nested dictionary, first level of key is the filenames, 
    # each filename has a subdict with 'hdr' and 'img' keys
    files = get_img_from_list(file_dir, filenames)
    return(files)



def search_imgs(metadata, variable, relation, value):
    '''PURPOSE: Look into metadata for images with the assigned characteristic.
        INPUT:
            metadata (dataframe) = Metadata of the group of images among which the search must be performed. 
            variable (string) = Variable to consider for selecting images. It is the selection criteria together with value and relation
            relation (string) = Relation between the variable and the images' characteristic. Valid values are ['major', 'equal', 'minor']
            value (sting or int/float) = value of the variable to select images
        OUTPUT:
            subset_metadata (dataframe) = Contain metadata of images with the assigned characteristic
        INVOKED FUNCTIONS:
            None'''
    # TO DO: check the variable is present
    if relation in ['major', 'equal', 'minor']: 
        if relation == 'equal':
            subset_metadata = metadata[metadata.loc[:,variable]==value]
        elif relation == 'major':
            subset_metadata = metadata[metadata.loc[:,variable]>value]
        elif relation == 'minor':
            subset_metadata = metadata[metadata.loc[:,variable]<value]
    else: 
        # Rise error if the assigned relation is not valid 
        raise ValueError('Assigned relation invalid')
    return(subset_metadata) 


def query_imgs(metadata, variable, relation, value):
    '''PURPOSE: Look into metadata for images with the assigned characteristic and return filenames.
        INPUT:
            metadata (dataframe) = Metadata of the group of images among which the search must be performed. 
            variable (string) = Variable to consider for selecting images. It is the selection criteria together with value and relation
            relation (string) = Relation between the variable and the images' characteristic. Valid values are ['major', 'equal', 'minor']
            value (sting or int/float) = value of the variable to select images
        OUTPUT:
            filenames (list) = List of filenames of the found images with the assigned characteristic
        INVOKED FUNCTIONS:
            None'''
    # Search for images with the sassigned characteristics  
    search_result = search_imgs(metadata, variable, relation, value)
    # Extract filename list for images acquired at the current value of exp_time
    filenames = list(search_result.loc[:,'FILENAME'])  
    return(search_result,filenames)


def avg_imgs(data):
    ''' PURPOSE: Generate an image which is the average of the images provided in input.
        INPUT:
            data (dict) = contains images and the related metadata. Structure is the following:
                data = { [filename1.fits] = { ['img']=[…], ['hdr']=[…] }, [filename2.fits] = {…} }
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


def imgs_stats(data, X_sigma_clip):
    ''' PURPOSE: Generate images of the average, the std, and the median of the images provided in input.
        INPUT:
            data (dict) = contains images and the related metadata. Structure is the following:
                data = { [filename1.fits] = { ['img']=[…], ['hdr']=[…] }, [filename2.fits] = {…}, ... }
        OUTPUT:
            stats (dict) = All statistics related to the priovided group of images. Structure is the following:
                stats = { ['avg'] (2D array) = [...] ; ['std'] (2D array) = [...] ; ['var'] (2D array) = [...] ; ['median'] (2D array) = [...]; 
                          ['mad'] (2D array) = [...] ; ['delta'] (2D array) = [...] ; ['avg_clip'] (2D array) = [...]} ; ['var_clip'] (2D array) = [...]
        INVOKED FUNCTIONS:
            None'''
    imgs={}; stats={}
    # Extract images from dictionary
    for filename in list(data.keys()):
        imgs[filename] = data[filename]['img']
    # Stack images contained in the data dictionary to compute the mean
    stacked_imgs = np.stack(list(imgs.values()), axis=2)
    # Average the received images and save the result in the output dictionary
    stats['avg'] = np.mean(stacked_imgs, axis=2)
    # Standard deviation for the averaged image
    stats['std'] = np.std(stacked_imgs, axis=2)
    # Variance for the averaged image
    stats['var'] = np.var(stacked_imgs, axis=2)
    # Median
    stats['median'] = np.median(stacked_imgs, axis=2)
    # Median Absolute Deviation
    N = len(stacked_imgs[0,0,:])
    stats['mad'] = np.sum(abs(stacked_imgs-np.median(stacked_imgs, axis=-1, keepdims=True)), axis=-1)/N
    # Calculate the difference between the highest and the lowest value for each pixel in the group of images
    max_values = np.amax(stacked_imgs, axis=2)      # 2D arrays with the maximum value for each pixel
    min_values = np.amin(stacked_imgs, axis=2)      # 2D arrays with the minimum value for each pixel
    delta = max_values-min_values                   # Biggest difference among all measurements from the same pixel
    stats['delta'] = delta
    # Calculate the clipped average by excluding data with a difference bigger than X*sigma_clip, X is an input
    stats['mad'][np.where(stats['mad']==0)]=1*10**-10      # Avoid to divide with zeros
    # Evaluate sigma clip
    sigma_clip = abs(stacked_imgs-np.median(stacked_imgs, axis=-1, keepdims=True))/stats['mad'][..., np.newaxis]
    # For every pixel it looks for values which difference is higher than X_sigma_clip*sigma_clip and change the value with a nan
    stacked_imgs_clip = np.where(sigma_clip>=X_sigma_clip, np.nan, stacked_imgs)
    # Compute the clipped average (nan values are not considered)
    stats['avg_clip'] = np.nanmean(stacked_imgs_clip, axis=2)
    # Compute the variance
    stats['var_clip'] = np.nanvar(stacked_imgs_clip, axis=2)
    return(stats)


def subtract_imgs(img1, img2):
    '''PURPOSE: Subtract img2 from img1.
        INPUT:
            img1 (2D array) = first image
            img2 (2D array) = second image to be subtracted from img1
        OUTPUT:
            new_img (2D array) = result of the image subtraction (img1-img2)
        INVOKED FUNCTIONS:
            None'''
    new_img = np.subtract(np.array(img1,dtype=np.float32),np.array(img2,dtype=np.float32))
    new_img[new_img < 0] = 0
    return(new_img)


def gen_datacube(dataset, img_type):
    ''' PURPOSE: Create a datacube with the group of images provided as input.
        INPUT:
            dataset (dict) = It containes the images to be used for generating a datacube
        OUTPUT: 
            datacube (3D array) = pile of images
        INVOKED FUNCTIONS:
            None'''
    imgs = {}
    # Extract images from dictionary
    for filename in list(dataset.keys()):
        imgs[filename]=dataset[filename][img_type]
    # Stack images contained in the data dictionary to compute the mean
    datacube = np.stack(list(imgs.values()), axis=2)
    return(datacube)


def generate_csv_metadata(dataset, path):
    '''Collect metadata in a Pandas dataframe and save them in a csv file in the assigned directory, if any.
        INPUT: 
            dataset (dict) = fits file list, for each filename, data and header. Expected structure is files { ['filename']  {['img], ['hdr]} }
            path (string) = path where the metadata will be saved
        OUTPUT: 
            metadata (dataframe) = table with the metadata of the images provided as input
        INVOKED FUNCTIONS:
            None'''
    # Initialize the dataframe by adding the metadata of the first file in the list as first row. 
    first_filename = list(dataset.keys())[0]  
    metadata = pd.DataFrame(columns=list(dataset[first_filename]['hdr'].keys()))
    # Loop on all the images contained in the dictionary and fill the dataframe
    for image in dataset.keys():
        hdr = dataset[image]['hdr']
        hdr_dict = {k: [dict(hdr)[k],] for k in dict(hdr)}
        hdr_df = pd.DataFrame(hdr_dict)
        metadata = pd.concat([metadata, hdr_df], ignore_index = True)
        if path:
            filename = os.path.join(path, 'metadata.csv')
            # Save extracted metadata as csv file in the assigned directory
            metadata.to_csv(filename)
    return(metadata)


def save_fits(data, filename, path):
    ''' PURPOSE: Save data in input as a fits file composed by image and header
        INPUT:
            data (dict) = data and metadata to be stored into a fits file. Structure is the following: 
                data = { ['img']=[…], ['hdr']=[…] }
            filename ()
            path (str) = path pointing the directory where the fits file has to be saved
        OUTPUT:
            None
        INVOKED FUNCTIONS:
            None'''
    hdr = fits.Header()
    keys = list(data['hdr'].keys())
    for key in keys[7:]:     # TO DO: generalize this 7
        hdr[key] = data['hdr'][key]
    hdr['FILENAME'] = filename
    hdu = fits.PrimaryHDU(data=data['img'], header=hdr) 
    hdu.writeto(os.path.join(path, hdr['FILENAME']), overwrite=True, checksum=True)
    return


def get_pol_imgs(data):
    ''' PURPOSE: Create 4 different images, each one for a specific polarization
        INPUT:
            data (dict) = data and metadata to be stored into a fits file. Structure is the following: data {['img']:[], ['hdr']:[]}
        OUTPUT:
            data (dict) = four polarized images
        INVOKED FUNCTIONS:
            None'''
    image = data['img']
    # Separate 4 images and save the result in the same dictionary provided as input
    data['pol_1'] = image[::2, ::2]
    data['pol_2'] = image[::2, 1::2]
    data['pol_3'] = image[1::2, ::2]
    data['pol_4'] = image[1::2, 1::2]
    return(data)


def gen_stats(data, mask, detector):
    ''' PURPOSE: Collect all statistics for the provided image.
        INPUT:
            data (dict) = Polarized images to compute statistics, whole image and their metadata. Structure is the following: 
                data = {['img']:[], ['hdr']:[], ['pol_1']:[], ['pol_2']:[], ['pol_3']:[], ['pol_4']:[]}
            mask (dict) = Info on mask to be applied on images
            detector (dict) = dimension on the detector. 
        OUTPUT:
            stats (dict) = Calculated statistics for each polarized image. Structure is the following:
                stats = {['pol_1']:{['max']:[], ['min']:[], ['mean']:[], ['std']:[]}, ['pol_2']:{..}, ['pol_3']:{...}, ['pol_4']:{...}, ['metadata']:{[mask]:{..}} }
        INVOKED FUNCTIONS:
            - apply_mask
            - get_stats '''
    # Generate statistics
    pols = ['pol_1', 'pol_2', 'pol_3', 'pol_4']
    stats = {}
    for pol in pols:
        img = data[pol]
        img_masked = apply_mask(img, mask, detector)
        stats[pol] = get_stats(img_masked)
    # Add information on the mask used to select the area for evaluate statistics
    stats['metadata'] = {}
    stats['metadata']['mask'] = mask 
    return (stats)
    

def get_stats(img):
    ''' PURPOSE: Collect statistics on the image in input.
        INPUT:
            img (2D array) = image to compute statistics on
        OUTPUT:
            stats (dict) = contains statistics computed on the image provided as input. Structure is the following:
                stats = {['max']:[], ['min']:[], ['mean']:[], ['std']:[], ['N']:[]}  
        INVOKED FUNCTIONS:
            None.'''
    stats = {}
    stats['min'] = np.round(np.nanmin(img),2)
    stats['max'] = np.round(np.nanmax(img),2)
    stats['mean'] = np.round(np.nanmean(img),2)
    stats['std'] = np.round(np.nanstd(img),2)
    stats['N'] = np.count_nonzero(np.isnan(img))
    return(stats)


def create_mask(mask, detector):
    ''' PURPOSE: Create a donut-shape mask to exclude pixels in the corner and behind the external occulter.
        INPUT:
            mask (dict) = Information to create the mask
            detector (dict) = information on detector dimensions
        OUTPUT:
            mask_img (2D array) = image composed of 1 and Nan values; Nan values correspond to portion of the image to not be considered, 
                1 correspond to the pixels within the area of interest.
        INVOKED FUNCTIONS:
            None '''
    # Generate a 2D arrays of zero values with the same shape of the detector
    mask_img = np.zeros((detector['npixls_y'], detector['npixls_x']))
    rows = np.arange(0, detector['npixls_y'], 1)
    columns = np.arange(0, detector['npixls_x'], 1)
    # Loop on row and columns to replace zeros with ones in the are between the two circumferences 
    for row in rows:
        for column in columns: 
            distance_in = np.sqrt((mask['y_center']-row)**2+(mask['x_center']-column)**2)
            distance_out = np.sqrt((mask['y_center']-row)**2+(mask['x_center_ext']-column)**2)
            if distance_in>mask['r_min'] and distance_out<mask['r_max']:
                mask_img[row, column] = 1
    # Replace zeros with nan. TO DO: verify if this is necessary
    mask_img[mask_img == 0] = np.nan
    return(mask_img)


def apply_mask(image, mask, detector):
    ''' PURPUSE: Generate and apply mask to data.
        INPUT:
            image (2D array) = image to be masked
            mask (dict) = information to create the mask
            detector (dict) = information on detector physical dimensions
        OUTPUT:
            masked_img (2D array) = image with mask applied
        INVOKED FUNCTIONS:
            - create_mask '''
    mask_img = create_mask(mask, detector)
    # Apply mask to image
    image_masked = image*mask_img        
    return(image_masked)


def get_demod(path):
    ''' PURPUSE: Read the fits file with the demodulatation tensor
        INPUT:
            - path (str) = pointing to the directory where the demodulation tensor images are stored 
        OUTPUT:
            - demod (dict) = demodulation tensor images. The structure is:
                demod = {['d_00'] = [...], ['d_01'] = [...], ['d_02'] = [...], ...}
        INVOKED FUNCTIONS: 
            - get_img_from_dir'''
    demod = {}
    # Read all fits files in the directory and save them as a dictionary 
    demod_files = get_img_from_dir(path)
    # Derive the element name
    filenames = list(demod_files.keys())
    demod = dict((k.removesuffix('.fits') , demod_files[k]['img']) for k in filenames)
    return(demod)


def rebin(arr, new_shape):
    ''' PURPUSE: Bin the assigned image 
        INPUT:
            - arr (2D array) = image to bin
            - new_shape (list) = shape of the output image, element [0] refers to rows, element [1] refers to columns
        OUTPUT:
            - binned (2D array) =  binned image 
        INVOKED FUNCTIONS:
            None'''
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    binned = arr.reshape(shape).mean(-1).mean(1)
    return(binned)


def remove_bad_pxls(data, bad_pxl_map, path=None):
    ''' PURPOSE: 
        INPUT:
        OUTPUT:
        INVOKED FUNCTIONS: '''
    # Initialize dictionary
    data_no_bad_pxl ={}
    # Generate metadata
    # TO DO: include metadata
    hdr = data['hdr']
    # Set bad pixel value as NaN (originally equal to -999.99)
    bad_pxl_map['img'][np.where(bad_pxl_map['img']==-999.999)] = np.nan
    data_no_bad_pxl['hdr'] = hdr 
    data_no_bad_pxl['img'] =  data['img']*bad_pxl_map['img']
    return(data_no_bad_pxl)


def get_tB(data):
    ''' PURPOSE: derive the image of total brightness by using the 4 polarized
                portions in the provided image.
        INPUT:
            - img (dict) = contains the image for the total brightness
                computation
        OUTPUT:
            - tB (2D array) = total brightness image
        INVOKED FUNCTIONS:
            - get_pol_imgs'''
    # Divide the full image into the 4 polarized images
    data = get_pol_imgs(data)
    # Evaluate the total brightness
    # TO DO: consider different pixels QE
    data['tB'] = (data['pol_1']+data['pol_2']+data['pol_3']+data['pol_4'])/2
    return(data)


def plot_img(image, title, min_val=None, max_val=None, path=None):
    '''Display an image acquired by CODEX.
    INPUT:
        - image (3D array) = 
        - title (str) = 
        - min_val = 
        - max_val  = 
        - path (str)= '''
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18}
    fig, ax = plt.subplots(figsize=(14, 10), sharey=True)
    fig.patch.set_facecolor('white')
    if min_val and max_val: 
        img=ax.imshow(image, vmin=min_val, vmax=max_val,  cmap='gray')
    else:
        img=ax.imshow(image, vmin=image.min(), vmax=image.max(),  cmap='gray')
    ax.set_title(title,  fontdict=font)
    ax.set_xlabel('Pixel #',  fontdict=font)
    ax.set_ylabel('Pixel #',  fontdict=font)
    fig.colorbar(img)
    if path:
        fig.savefig(os.path.join(path, title+'.png'))
    return


def plot_img_mask(image, title, mask, path=None, filename=None):
    '''Plot an image with mask boundaries'''
    fig, ax = plt.subplots(figsize=(14, 10), sharey=True)
    fig.patch.set_facecolor('white')
    fig.suptitle(title, fontsize=20, y=0.95)
    out_radius = plt.Circle((mask['x_center_ext'], mask['y_center']), mask['r_max'], color='red', linewidth=2.5, fill=False)
    in_radius = plt.Circle((mask['x_center'], mask['y_center']), mask['r_min'], color='red', linewidth=2.5, fill=False) 
    ax.imshow(image, vmin=image.min(), vmax=image.max(),  cmap='gray')
    ax.add_patch(out_radius)
    ax.add_patch(in_radius)
    if path:
        fig.savefig(os.path.join(path, filename+'.png'))
    return

def plot_pxl_distr(img, title, x_label, y_label):
    ''' PURPOSE: 
        INPUT:
            stats (dict) = Statistics computed over the image. Keys are pixel values and value is the number of pixels with that value 
            title (str) = Title for the plot
            x_label (str) = x-axis label
            y_label (str) = y-axis label
        OUTPUT:
        INVOKED FUNCTIONS: '''
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18}
    # Prepare dataset to be plotted as an histogram: from 2D to 1D array
    flat_data = img.flatten()
    fig, ax = plt.subplots(figsize=(14, 10), sharey=True)
    hist(flat_data, bins='freedman')
    ax.grid()
    ax.set_title(title,  size =20)
    ax.set_xlabel(x_label,  fontdict=font)
    ax.set_ylabel(y_label,  fontdict=font)
    return


def plot_img_area(image, title, x_ul, y_ul, width, heigh, path, filename):
    '''Plot an image with mask boundaries'''
    fig, ax = plt.subplots(figsize=(14, 10), sharey=True)
    fig.patch.set_facecolor('white')
    fig.suptitle(title, fontsize=20, y=0.95)
    rect = patches.Rectangle((x_ul, y_ul), width, heigh, linewidth=1, edgecolor='r', facecolor='none') 
    ax.imshow(image, vmin=image.min(), vmax=image.max(),  cmap='gray')
    ax.add_patch(rect)
    #fig.savefig(os.path.join(path, filename+'.png'))
    return


def plot_slicex(target, title, xslice, yslice_min, yslice_max, path, filename):
    '''Generate a plot with and image and a image profile for an assigned position'''
    # Plot style
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18}
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(14, 6)) 
    ax0.imshow(target, vmin=target.min(), vmax=target.max(), cmap='gray', aspect="auto")
    ax0.plot([xslice,xslice],[yslice_min,yslice_max], color='r', alpha=0.5)
    ax0.set_title('Vertical cut', fontdict=font)
    ax1.plot(target[yslice_min:yslice_max,xslice], '0.7', color='r')
    ax1.set_title('Profile',  fontdict=font)
    ax1.set_xlabel('Pixel #',  fontdict=font)
    ax1.set_ylabel('DN',  fontdict=font)
    ax1.grid()
    fig.suptitle(title, size=20, y=1.1)
    fig.tight_layout()
    fig.savefig(os.path.join(path, filename+'.png'))
    return


def plot_slicey(target, title, xslice_min, xslice_max, yslice, path, filename):
    '''Generate a plot with and image and a image profile for an assigned position'''
    # Plot style
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18}
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(14, 6)) 
    ax0.imshow(target, vmin=target.min(), vmax=target.max(), cmap='gray', aspect="auto")
    ax0.plot([xslice_min,xslice_max],[yslice,yslice], color='r', alpha=0.5)
    ax0.set_title('Horizontal cut', fontdict=font)
    ax1.plot(target[yslice,xslice_min:xslice_max], '0.7', color='r')
    ax1.set_title('Profile', fontdict=font)
    ax1.set_xlabel('Pixel #',  fontdict=font)
    ax1.set_ylabel('DN',  fontdict=font)
    ax1.grid()
    fig.suptitle(title, size=20, y=1.1)
    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    fig.savefig(os.path.join(path, filename+'.png'))
    return


def plot_pol_imgs(images, title):
    '''Plot the 4 images with different polarization'''
    # images = normalize_images(images)
    fig, axs = plt.subplots(2,2, figsize=(15, 11), sharey=True)
    fig.patch.set_facecolor('white')
    fig.suptitle(title, fontsize=20, y=0.95)
    img1 = axs[0,0].imshow(images['pol_1'], vmin=0, vmax=4096, cmap='Reds', aspect="auto")
    axs[0,0].set_title('Pol 1', fontsize=16)
    fig.colorbar(img1, ax=axs[0,0])
    img2 = axs[0,1].imshow(images['pol_2'], vmin=0, vmax=4096, cmap='Reds', aspect="auto")
    axs[0,1].set_title('Pol 2', fontsize=16)
    fig.colorbar(img2, ax=axs[0,1])
    img3 = axs[1,0].imshow(images['pol_3'], vmin=0, vmax=4096, cmap='Reds', aspect="auto")
    axs[1,0].set_title('Pol 3', fontsize=16)
    fig.colorbar(img3, ax=axs[1,0])
    img4 = axs[1,1].imshow(images['pol_4'], vmin=0, vmax=4096, cmap='Reds', aspect="auto")
    axs[1,1].set_title('Pol 4', fontsize=16)
    fig.colorbar(img4, ax=axs[1,1])
    return


def normalize_images(images):
    maxs=[]; mins=[]
    for image in images.keys(): 
        maxs.append(np.max(images[image]))
        mins.append(np.min(images[image]))
    Max = np.max(maxs)
    Min = np.min(mins)
    norm_images={}
    for image in images.keys(): 
        norm_images[image] = (images[image]-np.min(images[image]))*((Max-Min)/(np.max(images[image])-np.min(images[image])))+Min
    return(norm_images)


def compose_image(images):
    new_image = np.zeros((3000, 4096))
    for image_id in images.keys():
        ymin = images[image_id]['ymin']
        xmin = images[image_id]['xmin']
        ymax = images[image_id]['ymax']
        xmax = images[image_id]['xmax']
        new_image[ymin:ymax, xmin:xmax] = images[image_id]['image'][ymin:ymax, xmin:xmax]
    return(new_image)


def plot_data(x, x_label, y, y_label, title, path, filename):
    '''Generates plot for data serie'''
    # Plot style
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }
    # Create plot
    plt.figure(figsize=(10,6))
    plt.plot(x,y, 'o')
    plt.title(title, fontdict=font)
    plt.xlabel(x_label, fontdict=font)
    plt.ylabel(y_label, fontdict=font)
    plt.grid(True)
    plt.savefig(os.path.join(path, filename+'.png'))
    return


def double_plot_data(x, x_label, y1, y1_label, y1_limit, y1_ref_value, y2, y2_label, title, path, filename):
    '''Generates plot for data serie'''
    # Plot style
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }
    # Create plot
    fig, ax1 = plt.subplots(figsize=(10,6))
    ax2 = ax1.twinx()
    ax1.axhline(y=y1_ref_value, color='r', alpha=0.5, linestyle='--', label = "Full Well Capacity")
    ax1.plot(x, y1, '-o', color='r')
    ax1.set_xlabel(x_label, fontdict=font)
    ax1.set_ylabel(y1_label, color='r', size=16, fontname='serif')
    ax1.set_ylim([0, y1_limit])
    ax1.tick_params(axis="y", labelcolor='r')
    ax1.grid(axis="x")
    ax2.plot(x, y2, '-o', color='b')
    ax2.set_xlabel(x_label,  fontdict=font)
    ax2.set_ylabel(y2_label, color='b', size=16, fontname='serif')
    ax2.tick_params(axis="y", labelcolor='b')
    ax2.grid(axis="x")
    fig.suptitle(title, fontdict=font, size=18)
    ax1.legend(loc='upper left')  
    plt.savefig(os.path.join(path, filename+'.png'))
    return



def plot_vign_profile(x, x_lab, y, y_lab, title, path, filename, invert_xaxis=0):
    # Plot style
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }
    # Create plot
    plt.figure(figsize=(10,6))
    plt.plot(x,y, '-')
    xticks= np.arange(0, round(max(x),2), 0.2)
    # Invert x, if necessary
    if invert_xaxis==1:
        plt.gca().invert_xaxis()
        plt.xticks(xticks[::-1])
    plt.title('Vignetting - Vertical Cut', fontdict=font)
    plt.ylim(0, 1.0)
    plt.xlabel('Degree of field [deg]', fontdict=font)
    plt.ylabel('Fraction of unvignetted rays', fontdict=font)
    plt.grid(True)
    return

def plot_pxl_stat(image, title, DN_sat, ref_value):
    '''Plots an histogram showing the statistics related to input image's pixels'''
    # Plot style
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }
    # Create plot
    plt.figure(figsize=(16,10))
    image_1d = image.flatten()
    n = plt.hist(image_1d[~np.isnan(image_1d)], 400)
    # pxls_distr = Counter(np.round(image[~np.isnan(image)]))
    # pxls_distr=sorted(pxls_distr.items())
    # DN, frequency = zip(*pxls_distr)
    # plt.bar(DN,frequency, 10.0)
    mean = np.nanmean(image_1d)
    stdev = np.nanstd(image_1d)
    plt.axvline(mean, color='r', linestyle='dashed', linewidth=1.2)
    plt.text(mean, np.max(n[0])*0.9, 'Mean: {:.2f}'.format(mean)+'\nStDev: {:.2f}'.format(stdev), color='r',fontsize=13)
    plt.axvline(mean, color='r', linestyle='dashed', linewidth=1.2)
    plt.text(mean, np.max(n[0])*0.9, 'Mean: {:.2f}'.format(mean)+'\nStDev: {:.2f}'.format(stdev), color='r',fontsize=13)
    plt.title(title, fontdict=font, size=18)
    plt.xlabel('DN', fontdict=font)
    plt.ylabel('Number of Pixels', fontdict=font)
    plt.grid(True)
    return

def plot_contour(data, detector, title):
    """ PURPOSE: Plot images applying contour.
        INPUT: 
            data (dict) = 
            detector (dict) = 
        OUTPUT:
            None
        INVOKED FUNCTIONS:
            None"""
    # Prepare the grid for data
    x = np.arange(0, detector['npixls_x'], 1)
    y = np.arange(0, detector['npixls_y'], 1)
    X, Y = np.meshgrid(x, y)
    # Prepare figure
    fig, axs = plt.subplots(2,2, figsize=(17, 11), sharey=True)
    fig.patch.set_facecolor('white')
    fig.suptitle(title, fontsize=20, y=0.95)
    pols = ['pol_1', 'pol_2', 'pol_3', 'pol_4']
    for pol in pols: 
        Z = data[pol] 
        CS3 = ax2.contourf(X, Y, Z)
    return