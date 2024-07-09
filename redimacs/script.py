from pathlib import Path
from astropy.io import fits
import argparse


from .combine_calibrations import (make_main_bias_from_directory,
                                   make_main_flat_from_directory)
from .apply_calibrations import (load_spectrum_and_apply_calibrations,
                                 load_image_and_apply_calibrations)
from .binning import bin_spectrum
from .sky_modelling import subtract_median
from .plots import show_image
from .wcs import generate_simple_wavelength_solution, save_fits_with_wcs


def reduce_spectrum():
    parser = argparse.ArgumentParser(description='Reduce IMACS long slit spectroscopy data.')
    parser.add_argument('raw_path', type=str, help='Path to the directory containing raw FITS files.')
    parser.add_argument('dataset_number', type=int, help='Dataset number to process.')
    parser.add_argument('--spectral_binning', type=int, help='Bin along the spectral axis', default=2)
    parser.add_argument('--save_name', type=str, help='Save name', default='')
    parser.add_argument('--lambda_min', type=float,
                        help='For rough wavelength map with WCS: lambda of leftmost pixel', default=None)
    parser.add_argument('--lambda_max', type=float,
                        help='For rough wavelength map with WCS: lambda of rightmost pixel', default=None)
    parser.add_argument('--no_sky_subtraction', dest='sky_subtraction', action='store_false',
                        help='Disable sky subtraction (default: enabled)')

    parser.set_defaults(sky_subtraction=True)

    args = parser.parse_args()
    directory = Path(args.raw_path)
    dataset_number = int(args.dataset_number)
    print(f'Reducing exposures {dataset_number} in directory {directory}')

    # process calibrations
    # bias
    make_main_bias_from_directory(directory)

    # flats
    make_main_flat_from_directory(directory, flat_type='spec')

    # Load, subtract sky, bin and save spectrum
    spectrum, headers = load_spectrum_and_apply_calibrations(directory, dataset_number)

    # save name
    if not args.save_name:
        # create a name using the coordinates
        try:
            ra = ''.join(headers[0]['ra'].strip().split(':')[:2])
            dec = ''.join(headers[0]['dec'].strip().split(':')[:2])
            name = f'J{ra}{dec}'
        except Exception as e:
            print(f'Naming of final file, error: {e}H. Saving as UNKNOWN')
            name = 'UNKNOWN'
    else:
        name = args.save_name

    if args.sky_subtraction:
        spectrum = subtract_median(spectrum)
    if args.spectral_binning > 1:
        spectrum = bin_spectrum(spectrum, args.spectral_binning)

    # Save the binned spectrum as a FITS file with WCS
    if args.lambda_min is not None and args.lambda_max is not None:
        image_shape = spectrum.shape
        wcs = generate_simple_wavelength_solution(image_shape, args.lambda_min, args.lambda_max)
        save_fits_with_wcs(directory / f'{name}_reduced_spectrum_{dataset_number:0>04}.fits',
                           spectrum, wcs, headers[0])

    # Show and save the image
    show_image(spectrum, directory / f'{name}_reduced_spectrum_{dataset_number:0>04}.jpg')


def reduce_image():
    parser = argparse.ArgumentParser(description='Reduce IMACS imaging data.')
    parser.add_argument('raw_path', type=str, help='Path to the directory containing raw FITS files.')
    parser.add_argument('dataset_number', type=int, help='Dataset number to process.')
    parser.add_argument('--save_name', type=str, help='Save name', default='')
    parser.add_argument('--no_plate_solving', dest='plate_solving', action='store_false',
                        help='Disable plate solving (default: enabled)')
    parser.add_argument('--redo', dest='redo', action='store_false',
                        help='redo even if final file already exists? (default: do not redo)')
    parser.add_argument('--ccd_number', dest='ccd_number', type=int, default=None,
                        help='process only the given CCD number (1-8). Default: do all when not specified.')
    args = parser.parse_args()
    directory = Path(args.raw_path)
    dataset_number = int(args.dataset_number)
    print(f'Reducing exposures {dataset_number} in directory {directory}')

    # process calibrations
    # bias
    make_main_bias_from_directory(directory)

    # flats
    make_main_flat_from_directory(directory, flat_type='imag')

    # Load, subtract sky, bin and save spectrum
    ccds, headers = load_image_and_apply_calibrations(directory, dataset_number)

    # save name
    if not args.save_name:
        # create a name using the coordinates
        try:
            if type(headers) is dict:
                index = list(headers.keys())[0]
            else:
                index = 0
            ra = ''.join(headers[index]['ra'].strip().split(':')[:2])
            dec = ''.join(headers[index]['dec'].strip().split(':')[:2])
            name = f'J{ra}{dec}'
        except Exception as e:
            print(f'Naming of final file, error: {e}H. Saving as UNKNOWN')
            name = 'UNKNOWN'
    else:
        name = args.save_name

    # Save the images
    for key in ccds.keys():
        if args.ccd_number is not None and key != args.ccd_number:
            continue
        save_name = directory / f"{name}_c{key}_iff{args.dataset_number:0>4}.fits"
        if save_name.exists() and not args.redo:
            print(f'{save_name} already processed and no redo flag: skip.')
            continue
        array = ccds[key]
        header = headers[key]
        hdu = fits.PrimaryHDU(array, header)

        hdu.writeto(save_name, overwrite=True)
        if args.plate_solving:
            from widefield_plate_solver import plate_solve
            ra_approx = header['RA-D']
            dec_approx = header['DEC-D']
            binning = int(header['binning'].split('x')[0])
            wcs = plate_solve(save_name, sources=None, use_api=True, ra_approx=ra_approx, dec_approx=dec_approx,
                              scale_min=0.09, scale_max=0.12*binning, use_n_brightest_only=60)
            save_fits_with_wcs(save_name, array, wcs, header)

