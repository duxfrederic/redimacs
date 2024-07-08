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

    target_file = directory / f'iff{id_number:0>04}c{ccd_number}.fits'
    binning = fits.getheader(target_file).get('BINNING', '1x1')
    bx, by = binning.split('x')
    bx, by = int(bx), int(by)
    original_header = fits.getheader(target_file)
    slit = original_header.get('SLITMASK').lower()

    crop = 3600 // by
    crop2 = 100 // bx
    a = fits.getdata(target_file)
    a = a[crop:4096//by-crop2, :2048 // bx]
    b = fits.getdata(directory / f'stacked_bias_c{ccd_number}_{binning}.fits')
    b = b[crop:4096//by-crop2, :2048 // bx]

    f = fits.getdata(directory / f'stacked_flat_c{ccd_number}_{binning}_{slit}.fits')
    f = f[crop:4096//by-crop2, :2048 // bx]

    return (a - b) / f, bx, original_header


def load_spectrum_and_apply_calibrations(directory, id_number):
    ccds = []
    headers = []
    for ccd_number in [6, 5, 8, 7]:
        array, binning_x, original_header = load_ccd(directory, id_number, ccd_number)
        array = array[:, ::-1]
        ccds.append(array)
        gap_pixels = 57 // binning_x
        gap_array = np.zeros((array.shape[0], gap_pixels))
        ccds.append(gap_array)  # Add the gap array between CCDs
        headers.append(original_header)

    ccds.pop()  # Remove the last gap array as it's not needed after the last CCD

    # Concatenate all arrays along the short axis (axis 1)
    spectrum = np.concatenate(ccds, axis=1)

    # aaand inverse, because red is on the left on imacs
    return spectrum[:, ::-1], headers
