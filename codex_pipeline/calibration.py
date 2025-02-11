import os

from astropy.io import fits
import pandas as pd

def get_qe(path, wavelength):
    qe = {}
    pols = ['pol_1', 'pol_2', 'pol_3', 'pol_4']
    for pol in pols:
        qe[pol] = pd.read_csv(os.path.join(path, pol+'.csv'), index_col=[0])
    qe_wl = get_filter_qe(qe, wavelength, pols)
    return(qe_wl)


def get_filter_qe(qe, wavelength, pols):
    # Get quantum efficiency for this wavelength
    qe_wl = {}
    for pol in pols:
        qe_wl[pol] = qe[pol].loc[wavelength, 'mean']
    return(qe_wl)


def get_vign(path):
    ''' PURPOSE: Read and collect vignetting images in a dictionary.
        INPUT:
            - path (str) = points at the directory where the vignetting file is
                saved and stored
        OUTPUT:
            - vign (dict) = 4 images of vignetting for the 4 polarization imgs
        INVOKED FUNCTIONS:
            None '''
    # VIGNETTING
    # Vignetting filename. This filename might change to show
    # additional metadata as the version
    vign_filename = 'vignetting.fits'
    # Read and get content of vignetting fits file.
    # This fits file contains a vignetting function for each one of the four
    # polarized images. Data are within file extensions
    vignetting = fits.open(os.path.join(path, vign_filename),
                           checksum=True)
    vign = {}
    vign['pol_1'] = vignetting[1].data
    vign['pol_2'] = vignetting[2].data
    vign['pol_3'] = vignetting[3].data
    vign['pol_4'] = vignetting[4].data
    return(vign)
