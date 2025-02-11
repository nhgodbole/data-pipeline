"""
CODEX Polarization function.
Created on Mon Jan  6 19:14:14 2025

@author:
    mcasti
    The Catholic University of America at NASA/GSFC
    email: marta.casti@nasa.gov
"""

import os
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable

from codex_pipeline import *


# __all__ = ["get_avgd_img", "analyze_pol", "derive_stokes", "get_aolp",
#            "get_dolp", "plot_stokes", "plot_aolp", "plot_dolp", "save_stokes_fits"]


# def get_avgd_img(data, metadata, feat):
#     ''' PURPOSE: Compute the averaged image of the images provided as input.
#         INPUT:
#             data (dict) = contains images and the related metadata 
#             metadata (dataframe) = metadata related to the images in data. These are the UPDATED metadata
#             feat (string) = to be reported in the filename of the averaged image  
#         OUTPUT:
#             avg_img (dict) = contains the image resulting from the average and its metadata. Structure is avg_img{ [filename] {['img'], ['hdr']} }
#         INVOKED FUNCTIONS:
#             - avg_imgs '''
#     # Initialized output dictionary 
#     avg_img = {}
#     # IMAGE
#     # Evaluate the average of the images and the associated standard deviation
#     avg_img['img'], avg_img['std'] = avg_imgs(data)
#     # HEADER
#     # Generate new header for the resulting image starting from the metadata of the first considered image
#     new_header = metadata.iloc[0].to_dict()
#     # Remove first key of the dictionary, since it is the pandas dataframe index (TO DO: Verify)
#     first_key=list(new_header.keys())[0]
#     del(new_header[first_key])
#     # Generate a filename for the new image
#     new_header['FILENAME'] = feat+'_avg'
#     # Create 'HISTORY' keyword
#     new_header['HISTORY'] = 'Result of the average of images'
#     # TO DO: add the acquisition time of the first and the last images used to create the averaged image
#     #for n in np.arange(1,6,1):
#     #    new_header['TEMP'+str(n)]=np.mean(metadata.loc[:,'TEMP'+str(n)])
#     # Save new header in the dictionary
#     avg_img['hdr'] = new_header
#     return(avg_img)

def get_pol_imgs(data):
    ''' PURPOSE: Separate the 4 polarized images
        INPUT:
            - data (dict) = data and metadata to be stored into a fits file.
                Structure is the following: data {['img']:[], ['hdr']:[]}
        OUTPUT:
            data (dict) = four polarized images
        INVOKED FUNCTIONS:
            None'''
    image = data['img']
    # Separate 4 images and save the result in the input dictionary
    data['pol_1'] = image[::2, ::2]
    data['pol_2'] = image[::2, 1::2]
    data['pol_3'] = image[1::2, ::2]
    data['pol_4'] = image[1::2, 1::2]
    return(data)


def get_demod(demod_path):
    ''' PURPOSE: Read fits files containing the demodulation tensor images and
            save content in a dictionary
        INPUT:
            - demod_path (string) = Directory that contains the demodulation
                tensor images
        OUTPUT:
            - demod (dict) = collects the images of the demodulation tensor
                elements
        INVOKED FUNCTIONS:
            - util.get_img_from_dir'''
    # Initialize dictionary to contain the demodulation tensor elements
    demod = {}
    # Read all fits files in the directory and save them as a dictionary
    demod_files = get_img_from_dir(demod_path)
    # Derive the element name
    filenames = list(demod_files.keys())
    demod = dict((k.removesuffix('.fits'), demod_files[k]['img']) for k in
                 filenames)
    return(demod)


def analyze_pol(img, demod, detector, mask, gain, qe, plot=0, path=None):
    ''' PURPOSE: Execute the polarization analysis of the input image.
        INPUT:
            - img (dict) = contains the image to be analyzed and its header
            - demod (dict) = contains the images of the demodulation tensor
            - detector (dict) = info related on the size of the detector
            - gain (float) = CODEX detector gain
            - qe (dic) = filter measured quantum efficiency
            - plot (bool) = flag to plot or not
            - path (string) =
        OUTPUT:
            - stokes (dict) = Pixel-by-pixel map with the derived Stokes vector
                elements
            - aolp (array 2D) = Pixel-by-pixel map with the derived Angle of
                linear polarization
            - dolp (array 2D) = Pixel-by-pixel map with the derived Degree of
                linear polarization
        INVOKED FUNCTIONS:
            - polarization.derive_stokes
            - polarization.plot_stokes
            - polarization.plot_aolp
            - polarization.plot_dolp'''
    # TO DO: modify to integrate with the pipeline:
    #   -> add strategy to manage metadata
    #   ->
    # Derive Stokes vector
    stokes = derive_stokes(img, demod, detector, gain, qe)
    # Derive Angle of linear polarization
    aolp = 1/2*np.arctan2(stokes['S_2'], stokes['S_1'])
    # Derive degree of linear polarization
    dolp = np.sqrt((stokes['S_1'])**2+(stokes['S_2'])**2)/(stokes['S_0'])
    if plot == 1:
        # Generate plots and save them as images and as fits files
        plot_stokes(stokes, 'Stokes Vector', mask, detector, path,
                    'Stokes_Vector')
        plot_aolp(aolp, 'Angle of Linear Polarization ', mask, detector, path,
                  'AoLP')
        plot_dolp(dolp, 'Degree of Linear Polarization ', mask, detector, path,
                  'DoLP')
    if path:
        # Save fits file with Stokes vector
        hdr = fits.Header()
        hdr['STOKES_E'] = 'S_0'
        hdu = fits.PrimaryHDU(data=stokes['S_0'], header=hdr)
        hdu1 = fits.ImageHDU(stokes['S_1'])
        hdu1.header['STOKES_E'] = 'S_1'
        hdu2 = fits.ImageHDU(stokes['S_2'])
        hdu2.header['STOKES_E'] = 'S_2'
        new_hdul = fits.HDUList([hdu, hdu1, hdu2])
        new_hdul.writeto(os.path.join(path, 'StokesVector.fits'),
                         overwrite=True)
        # Save fits file with the Angle of Linear Polarization
        hdr = fits.Header()
        hdu = fits.PrimaryHDU(data=aolp, header=hdr)
        hdu.writeto(os.path.join(path, 'AoLP.fits'), overwrite=True,
                    checksum=True)
        # Save fits file with the Degree of Linear polarization
        hdr = fits.Header()
        hdu = fits.PrimaryHDU(data=dolp, header=hdr)
        hdu.writeto(os.path.join(path, 'DoLP.fits'), overwrite=True,
                    checksum=True)
    return(stokes, aolp, dolp)


def derive_stokes(img, demod, detector, gain, qe):
    ''' PURPOSE: To derive the Stokes vector associated to the observed scene.
        INPUT:
            - img (dict) = image to demodulate
            - demod (dict) = demodulation tensor
            - detector (dict) = detector dimensions
        OUTPUT:
            - stokes (dict) = images of the Stokes's vector elements
        INVOKED FUNCTIONS:
            - get_pol_imgs '''
    pols = ['pol_1', 'pol_2', 'pol_3', 'pol_4']
    # Divide the image in the 4 polarizad images
    # img = get_pol_imgs(img)
    # Initialize dictionary that will contain Stokes vector images
    stokes_elements = ['S_0', 'S_1', 'S_2']
    # Generate a dictionary that contains the images of Stokes vector elements
    stokes = {}
    stokes[stokes_elements[0]] = np.zeros((detector['npixls_y'],
                                           detector['npixls_x']))
    stokes[stokes_elements[1]] = np.zeros((detector['npixls_y'],
                                           detector['npixls_x']))
    stokes[stokes_elements[2]] = np.zeros((detector['npixls_y'],
                                           detector['npixls_x']))
    # Normalize input data
    norm_data = {}
    for pol in pols:
        div = 0.5*(img['pol_1']/((1/gain)*qe['pol_1']) +
                   img['pol_2']/((1/gain)*qe['pol_2']) +
                   img['pol_3']/((1/gain)*qe['pol_3']) +
                   img['pol_4']/((1/gain)*qe['pol_4']))
        div[np.where(div == 0)] = 1 * 10 ** -10
        norm_data[pol] = (img[pol]/((1/gain)*qe[pol]))/div
    # Compute Stokes vector for each pixel within the image
    rows = np.arange(0, detector['npixls_y'], 1)
    columns = np.arange(0, detector['npixls_x'], 1)
    for row in rows:
        for column in columns:
            # For each pxl
            dem_matrix = np.array([[demod['d_00'][row, column],
                                    demod['d_01'][row, column],
                                    demod['d_02'][row, column],
                                    demod['d_03'][row, column]],
                                   [demod['d_10'][row, column],
                                    demod['d_11'][row, column],
                                    demod['d_12'][row, column],
                                    demod['d_13'][row, column]],
                                   [demod['d_20'][row, column],
                                    demod['d_21'][row, column],
                                    demod['d_22'][row, column],
                                    demod['d_23'][row, column]]])
            meas = np.array([norm_data['pol_1'][row, column],
                            norm_data['pol_2'][row, column],
                            norm_data['pol_3'][row, column],
                            norm_data['pol_4'][row, column]])
            meas = np.transpose(meas)
            stokes_par = dem_matrix.dot(meas)
            stokes['S_0'][row, column] = stokes_par[0]
            stokes['S_1'][row, column] = stokes_par[1]
            stokes['S_2'][row, column] = stokes_par[2]
    return(stokes)


def get_aolp(stokes):
    # Compute the Angle Of Linear Polarization of the image in input
    theta = 1/2*np.arctan2(stokes['S_2'], stokes['S_1'])
    return(theta)


def get_dolp(stokes):
    ''' PURPOSE: Generate theimage showing the degree of linear polarization
            of the observed scene.
        INPUT:
            - stokes (dict) = Images of the Stokes vector elements
        OUTPUT:
            - dolp (2D array) = Image of degree of linear polarization
        INVOKED FUNCTIONS:
            None'''
    # Substitute zeros with 1*10^-10 to avoid mathematical errors
    stokes['S_0'][np.where(stokes['S_0'] == 0)] = 1 * 10 ** -10
    # Derive the Degree of linear polarization
    dolp = np.sqrt((stokes['S_1'])**2+(stokes['S_2'])**2)/(stokes['S_0'])
    return(dolp)


def get_pb(stokes):
    ''' PURPOSE: Derive the image of Polarization brightness.
        INPUT:
            - stokes (dict) = Images of the Stokes vector elements
        OUTPUT:
            - pb (2D array) = Image of polarized brightness
            '''
    # Evaluate pB
    pb = np.sqrt((stokes['S_1']**2+stokes['S_2']**2))
    return (pb)


def plot_stokes(stokes, title, mask, detector, path=None, filename=None):
    ''' PURPOSE: Plot the derived Stokes vector for each pixel in the image
        INPUT:
            stokes (dict) =
            title (str) = Plot title
            mask (dict) = Mask parameters such as radii and centers
            detector (dict) = Detector parameters
            path (str) = Optional. Path pointing the folder where the plot will
                be save as png image
            filename (str) = Optional
        OUTPUT:
        INVOKED FUNCTIONS: '''
    # Margin from theo
    # theo_space = 0.1
    # Evaluate statistics - S0
    masked_stokes_im = apply_mask(stokes['S_0'], mask, detector)
    S0_distr = Counter(np.round(masked_stokes_im[~np.isnan(masked_stokes_im)],
                                1))
    mean = masked_stokes_im[~np.isnan(masked_stokes_im)].mean()
    stdev = masked_stokes_im[~np.isnan(masked_stokes_im)].std()
    # Start generating the plot
    fig, axs = plt.subplots(3, 2, figsize=(17, 18))
    txt_space = 0.02
    kwargs = {'format': '%.1f'}
    plt.suptitle(x=0.5, y=0.95, t=title, fontsize=26, ha="center",
                 transform=fig.transFigure)
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    fig.patch.set_facecolor('white')
    # Plot S0 image
    im1 = axs[0, 0].imshow(stokes['S_0'], vmax=mean+3*stdev, vmin=mean-3*stdev,
                           cmap='Purples_r', aspect="auto")
    axs[0, 0].axes.get_xaxis().set_ticks([])
    axs[0, 0].axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(axs[0, 0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im1, cax=cax, **kwargs)
    cbar.set_ticks([mean-0.1, mean, mean+0.1])
    im1 = axs[0, 0].set_title('$S_0$', fontsize=18, pad=20)
    # Plot S0 stats considering the masked image
    # axs[0, 1].bar(S0_distr.keys(), S0_distr.values(), width = 0.01,
    #               color='blue')
    xbins = np.arange(mean-0.1, mean+0.1, 0.02)
    # xbins = np.arange(-(mean-(stdev*3)), (mean+(stdev*3)), 0.02)
    axs[0, 1].hist(masked_stokes_im[~np.isnan(masked_stokes_im)], bins=xbins,
                   histtype='step', fill=True, color="purple", label='Image')
    # Derive statistics , i.e., mean and stdev
    ymax = axs[0, 1].get_ylim()[1]
    axs[0, 1].axvline(mean, color='k', linestyle='dashed', linewidth=1.2)
    axs[0, 1].text(mean+txt_space, ymax*0.8, 'Mean: {:.2f}'.format(mean) +
                   '\nStDev: {:.2f}'.format(stdev), color='k', fontsize=11)
    axs[0, 1].set_xlim([mean-0.5, mean+0.5])
    axs[0, 1].set_xlabel("$S_0$ Value", fontsize=13)
    axs[0, 1].set_ylabel("Counts", fontsize=13)
    axs[0, 1].set_title('Values distribution', fontsize=18, pad=20)
    # Evaluate statistics - S1
    masked_stokes_im = apply_mask(stokes['S_1'], mask, detector)
    S1_distr = Counter(np.round(masked_stokes_im[~np.isnan(masked_stokes_im)],
                                2))
    mean = masked_stokes_im[~np.isnan(masked_stokes_im)].mean()
    stdev = masked_stokes_im[~np.isnan(masked_stokes_im)].std()
    # Plot S1 image
    im2 = axs[1, 0].imshow(stokes['S_1'], vmax=mean+5*stdev, vmin=mean-5*stdev,
                           cmap='Purples_r', aspect="auto")
    axs[1, 0].axes.get_xaxis().set_ticks([])
    axs[1, 0].axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(axs[1, 0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im2, cax=cax, **kwargs)
    im2 = axs[1, 0].set_title('$S_1$', fontsize=18, pad=20)
    # Plot S1 stats
    # axs[1,1].bar(S1_distr.keys(), S1_distr.values(), width = 0.01, color='blue')
    xbins = np.arange(mean-3*stdev, mean+3*stdev, 0.02)
    axs[1, 1].hist(masked_stokes_im[~np.isnan(masked_stokes_im)], bins=xbins,
                   histtype='step', fill=True, color="purple", label='Image')
    ymax = axs[1, 1].get_ylim()[1]
    axs[1, 1].axvline(mean, color='k', linestyle='dashed', linewidth=1.2)
    axs[1, 1].text(mean+txt_space, ymax*0.8, 'Mean: {:.2f}'.format(mean) +
                   '\nStDev: {:.2f}'.format(stdev), color='k', fontsize=11)
    axs[1, 1].set_xlim([mean-0.2, mean+0.2])
    axs[1, 1].set_xlabel("$S_1$ Value", fontsize=13)
    axs[1, 1].set_ylabel("Counts", fontsize=13)
    axs[1, 1].set_title('Values distribution', fontsize=18, pad=20)
    # Evaluate statistics - S2
    masked_stokes_im = apply_mask(stokes['S_2'], mask, detector)
    S2_distr = Counter(np.round(masked_stokes_im[~np.isnan(masked_stokes_im)],
                                2))
    mean = masked_stokes_im[~np.isnan(masked_stokes_im)].mean()
    stdev = masked_stokes_im[~np.isnan(masked_stokes_im)].std()
    # Plot S2 image
    im3 = axs[2, 0].imshow(stokes['S_2'], vmax=mean+5*stdev, vmin=mean-5*stdev,
                           cmap='Purples_r', aspect="auto")
    axs[2, 0].axes.get_xaxis().set_ticks([])
    axs[2, 0].axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(axs[2, 0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im3, cax=cax, **kwargs)
    im3 = axs[2, 0].set_title('$S_2$', fontsize=18, pad=20)
    # Plot S2 stats
    # axs[2,1].bar(S2_distr.keys(), S2_distr.values(), width = 0.06, color='blue')
    xbins = np.arange(mean-3*stdev, mean+3*stdev, 0.02)
    axs[2, 1].hist(masked_stokes_im[~np.isnan(masked_stokes_im)], bins=xbins,
                   histtype='step', fill=True, color="purple", label='Image')
    ymax = axs[2, 1].get_ylim()[1]
    axs[2, 1].axvline(mean, color='k', linestyle='dashed', linewidth=1.2)
    axs[2, 1].text(mean+txt_space, ymax*0.8, 'Mean: {:.2f}'.format(mean) +
                   '\nStDev: {:.2f}'.format(stdev), color='k', fontsize=11)
    axs[2, 1].set_xlim([mean-0.2, mean+0.2])
    axs[2, 1].set_xlabel("$S_2$ Value", fontsize=13)
    axs[2, 1].set_ylabel("Counts", fontsize=13)
    axs[2, 1].set_title('Values Distribution', fontsize=18, pad=20)
    if path:
        plt.savefig(os.path.join(path, filename + '.png'))
        plt.savefig(os.path.join(path, filename + '.eps'))
    return


def plot_aolp(aolp, title, mask, detector, path=None, filename=None):
    ''' PURPOSE:
        INPUT:
        OUTPUT:
        INVOKED FUNCTIONS: '''
    # Evaluate statistics in the image
    masked_aolp = apply_mask(np.degrees(aolp), mask, detector)
    aolp_distr = Counter(np.round(masked_aolp[~np.isnan(masked_aolp)], 1))
    mean = masked_aolp[~np.isnan(masked_aolp)].mean()
    stdev = masked_aolp[~np.isnan(masked_aolp)].std()
    # Generate Plot
    fig, axs = plt.subplots(1, 2, figsize=(17, 5))
    fig.subplots_adjust(wspace=0.3)
    fig.patch.set_facecolor('white')
    plt.suptitle(x=0.5, y=1.1, t=title, fontsize=24, ha="center",
                 transform=fig.transFigure)
    im1 = axs[0].imshow(np.degrees(aolp), vmin=mean-(stdev*3),
                        vmax=mean+(stdev*3), cmap='Purples_r', aspect="auto")
    axs[0].axes.get_xaxis().set_ticks([])
    axs[0].axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.07)
    plt.colorbar(im1, cax=cax)
    im1 = axs[0].set_title('Angle of Linear Polarization', fontsize=18, pad=20)
    # Evaluate statistics in the image and plot density
    axs[1].bar(aolp_distr.keys(), aolp_distr.values(), width=0.1, color='blue')
    axs[1].axvline(mean, color='k', linestyle='dashed', linewidth=1.2)
    ymax = axs[1].get_ylim()[1]
    axs[1].text(mean+(stdev/2), ymax*0.9, 'Mean: {:.2f}'.format(mean) +
                '\nStDev: {:.2f}'.format(stdev), color='r', fontsize=11)
    axs[1].set_xlim([mean-(4*stdev), mean+(4*stdev)])
    axs[1].set_xlabel("AoLP Value", fontsize=15)
    axs[1].set_ylabel("Counts", fontsize=15)
    axs[1].set_title('Values Distribution', fontsize=18, pad=20)
    plt.subplots_adjust(top=0.8)
    if path:
        plt.savefig(os.path.join(path, filename+'.png'))
        plt.savefig(os.path.join(path, filename+'.eps'))
    return


def plot_dolp(dolp, title, mask, detector, path=None, filename=None):
    ''' PURPOSE:
        INPUT:
        OUTPUT:
        INVOKED FUNCTIONS: '''
    # Evaluate statistics in the image 
    masked_dolp = apply_mask(dolp, mask, detector)
    dolp_distr = Counter(np.round(masked_dolp[~np.isnan(masked_dolp)], 2))
    mean = masked_dolp[~np.isnan(masked_dolp)].mean()
    stdev = masked_dolp[~np.isnan(masked_dolp)].std()
    # Generate plot
    fig, axs = plt.subplots(1, 2, figsize=(17, 5))
    fig.subplots_adjust(wspace=0.3)
    plt.suptitle(x=0.5, y=1.1, t=title, fontsize=24, ha="center",
                 transform=fig.transFigure)
    im1 = axs[0].imshow(dolp, vmin=mean-(stdev*3), vmax=mean+(stdev*3),
                        cmap='jet', aspect="auto")
    axs[0].axes.get_xaxis().set_ticks([])
    axs[0].axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.07)
    plt.colorbar(im1, cax=cax)
    im1 = axs[0].set_title('Degree of Linear Polarization', fontsize=18,
                           pad=20)
    # Plot density distribution
    xbins = np.arange(-(mean-(stdev*3)), (mean+(stdev*3)), 0.02)
    axs[1].hist(masked_dolp[~np.isnan(masked_dolp)], bins=xbins, density=False,
                histtype='step', fill=True, color="purple", label='Image')
    axs[1].axvline(mean, color='k', linestyle='dashed', linewidth=1.2)
    ymax = axs[1].get_ylim()[1]
    axs[1].text(mean+(stdev/2), ymax*0.8, 'Mean: {:.2f}'.format(mean) +
                '\nStDev: {:.2f}'.format(stdev), color='k', fontsize=11)
    axs[1].set_xlim([0, 1.3])
    axs[1].set_xlabel("DoLP Value", fontsize=13)
    axs[1].set_ylabel("Density", fontsize=13)
    axs[1].set_title('Density distribution', fontsize=18, pad=20)
    if path:
        plt.savefig(os.path.join(path, filename + '.png'))
        plt.savefig(os.path.join(path, filename + '.eps'))
    return


def plot_pB(pB):
    fig, ax = plt.subplots(figsize=(15, 8), sharey=True)
    fig.patch.set_facecolor('white')
    img = ax.imshow(np.flipud(pB), vmin=0, vmax=1, cmap='Purples_r',
                    aspect='auto')
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    fig.colorbar(img)


def save_stokes_fits(stokes, path):
    ''' PURPOSE: Save the derived Stokes vector as a fits file.
        INPUT: 
            - stokes (dict): Contains the pixel-by-pixel derived Stokes vector elements.
        OUTPUT: 
            - 
        INVOKED FUNCTIONS:
            None'''
    hdr = fits.Header()
    hdr['STOKES_ELM'] = 'S_0'
    hdu = fits.PrimaryHDU(data=stokes['S_0'], header=hdr)
    hdu1 = fits.ImageHDU(stokes['S_1'])
    hdu1.header['STOKES_ELM'] = 'S_1'
    hdu2 = fits.ImageHDU(stokes['S_2'])
    hdu2.header['STOKES_ELM'] = 'S_2'
    new_hdul = fits.HDUList([hdu, hdu1, hdu2])
    new_hdul.writeto(os.path.join(path, 'StokesVector.fits'), overwrite=True)
    return