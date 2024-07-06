from astropy.wcs import WCS
from astropy.io import fits


def generate_wcs(image_shape, lambda_start=4000, lambda_end=8500):
    """
    Generate a WCS object for a 2D image with spatial (y) and spectral (x) axes.

    Parameters:
    - image_shape: tuple
        Shape of the 2D image (num_rows, num_cols).
    - lambda_start: float
        Wavelength corresponding to the first pixel in the spectral (x) direction.
    - lambda_end: float
        Wavelength corresponding to the last pixel in the spectral (x) direction.

    Returns:
    - wcs: WCS object
        The WCS object with the desired configuration.
    """
    num_rows, num_cols = image_shape

    w = WCS(naxis=2)
    
    # Set the reference pixel to the center of the image
    w.wcs.crpix = [num_cols / 2, num_rows / 2]
    
    # Set the reference coordinates
    w.wcs.crval = [(lambda_start + lambda_end) / 2, num_rows / 2]
    w.wcs.cdelt = [(lambda_end - lambda_start) / (num_cols - 1), 1]
    w.wcs.ctype = ['LINEAR', 'PIXEL']
    
    return w


def save_fits_with_wcs(filename, data, wcs):
    hdu = fits.PrimaryHDU(data, header=wcs.to_header())
    hdu.writeto(filename, overwrite=True)


