from astropy.io import fits
import numpy as np


def load_ccd(directory, id_number, ccd_number):
    """

    Args:
        directory:
        id_number:
        ccd_number:

    Returns:
        pre-red single CCD
    """
    crop = 3600
    crop2 = 100
    target_file = directory / f'iff{id_number:0>04}c{ccd_number}.fits'
    binning = fits.getheader(target_file).get('BINNING', '1x1')
    slit = fits.getheader(target_file).get('SLITMASK').lower()

    a = fits.getdata(target_file)
    a = a[crop:4096-crop2, :2048]
    b = fits.getdata(directory / f'stacked_bias_c{ccd_number}_{binning}.fits')
    b = b[crop:4096-crop2, :2048]

    f = fits.getdata(directory / f'stacked_flat_c{ccd_number}_{binning}_{slit}.fits')
    f = f[crop:4096-crop2, :2048]

    return (a - b) / f


def load_spectrum_and_apply_calibrations(directory, id_number):
    ccds = []
    gap_pixels = 92

    for ccd_number in [6, 5, 8, 7]:
        array = load_ccd(directory, id_number, ccd_number)[:, ::-1]
        ccds.append(array)
        gap_array = np.zeros((array.shape[0], gap_pixels))
        ccds.append(gap_array)  # Add the gap array between CCDs

    ccds.pop()  # Remove the last gap array as it's not needed after the last CCD

    # Concatenate all arrays along the short axis (axis 1)
    spectrum = np.concatenate(ccds, axis=1)
    return spectrum
