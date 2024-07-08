from pathlib import Path
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import argparse


from .combine_calibrations import (make_main_bias_from_directory,
                                   make_main_flat_from_directory)
from .apply_calibrations import load_spectrum_and_apply_calibrations
from .binning import bin_spectrum
from .sky_modelling import subtract_median
from .plots import show_image
from .wcs import generate_wcs, save_fits_with_wcs


def main():
    parser = argparse.ArgumentParser(description='Reduce IMACS long slit spectroscopy data.')
    parser.add_argument('raw_path', type=str, help='Path to the directory containing raw FITS files.')
    parser.add_argument('dataset_number', type=int, help='Dataset number to process.')
    parser.add_argument('--spectral_binning', type=int, help='Bin along the spectral axis', default=2)

    parser.add_argument('--no-sky_subtraction', dest='sky_subtraction', action='store_false',
                        help='Disable sky subtraction (default: enabled)')

    parser.set_defaults(sky_subtraction=True, flat_correction=True, bias_correction=True)

    args = parser.parse_args()
    directory = Path(args.raw_path)
    dataset_number = int(args.dataset_number)
    print(f'Reducing exposures {dataset_number} in directory {directory}')

    # process calibrations
    # bias
    make_main_bias_from_directory(directory)

    # flats
    make_main_flat_from_directory(directory)

    # Load, subtract sky, bin and save spectrum
    spectrum, headers = load_spectrum_and_apply_calibrations(directory, dataset_number)

    # create a name using the coordinates
    ra = ''.join(headers[0]['ra'].strip().split(':')[:2])
    dec = ''.join(headers[0]['dec'].strip().split(':')[:2])
    name = f'J{ra}{dec}'

    if args.sky_subtraction:
        spectrum = subtract_median(spectrum)
    if args.spectral_binning > 1:
        spectrum = bin_spectrum(spectrum, args.spectral_binning)

    # Save the binned spectrum as a FITS file with WCS
    image_shape = spectrum.shape
    wcs = generate_wcs(image_shape)
    save_fits_with_wcs(directory / f'{name}_reduced_spectrum_{dataset_number:0>04}.fits', spectrum, wcs, headers[0])

    # Show and save the image
    show_image(spectrum, directory / f'{name}_reduced_spectrum_{dataset_number:0>04}.jpg')


if __name__ == "__main__":
    main()

