from astropy.io import fits
import numpy as np


def load_ccd(directory, id_number, ccd_number, crop_bottom=0, crop_top=0):
    """

    Args:
        directory:
        id_number:
        ccd_number:
        crop_bottom: for spectra, object usually is in the middle, so we crop
        crop_top: same as above

    Returns:
        pre-red single CCD
    """

    target_file = directory / f'iff{id_number:0>04}c{ccd_number}.fits'
    binning = fits.getheader(target_file).get('BINNING', '1x1')
    bx, by = binning.split('x')
    bx, by = int(bx), int(by)
    original_header = fits.getheader(target_file)
    slit = original_header.get('SLITMASK').lower()


    a = fits.getdata(target_file)
    b = fits.getdata(directory / f'stacked_bias_c{ccd_number}_{binning}.fits')
    f = fits.getdata(directory / f'stacked_flat_c{ccd_number}_{binning}_{slit}.fits')

    if crop_bottom > 0 and crop_top > 0:
        crop = crop_bottom // by
        crop2 = crop_top // bx
        a = a[crop:4096//by-crop2, :2048 // bx]
        b = b[crop:4096//by-crop2, :2048 // bx]
        f = f[crop:4096//by-crop2, :2048 // bx]

    return (a - b) / f, (bx, by), original_header


def load_spectrum_and_apply_calibrations(directory, id_number):
    ccds = []
    headers = []
    for ccd_number in [6, 5, 8, 7]:
        array, binning, original_header = load_ccd(directory, id_number, ccd_number, crop_bottom=3700, crop_top=200)
        array = array[:, ::-1]
        ccds.append(array)
        gap_pixels = 167 // binning[0]
        print(f'Inserting gap: {gap_pixels}')
        gap_array = np.zeros((array.shape[0], gap_pixels))
        ccds.append(gap_array)  # Add the gap array between CCDs
        headers.append(original_header)

    ccds.pop()  # Remove the last gap array as it's not needed after the last CCD

    # Concatenate all arrays along the short axis (axis 1)
    spectrum = np.concatenate(ccds, axis=1)

    # aaand inverse, because red is on the left on imacs
    return spectrum[:, ::-1], headers


def load_image_and_apply_calibrations(directory, id_number):
    ccds = {}
    headers = {}
    for ccd_number in [1, 2, 3, 4, 6, 5, 8, 7]:
        array, binning, original_header = load_ccd(directory, id_number, ccd_number, crop_bottom=3700, crop_top=200)
        array = array[:, ::-1]
        ccds[ccd_number] = array
        headers[ccd_number] = original_header

    return ccds, headers
