import ccdproc
from astropy.nddata import CCDData
from astropy.io import fits
from pathlib import Path
import numpy as np

from .files_handling import get_list_of_files_in_directory


def stack_bias(directory, ccd_id, binning, read_mode, files, redo=False):
    combined_bias_file = directory / f'stacked_bias_{ccd_id}_{binning}_{read_mode}.fits'
    if combined_bias_file.exists() and not redo:
        return
    print(f'Combining {len(files)} for ccd {ccd_id}, binning {binning} and read mode {read_mode}')
    ccd_list = []
    for file in files:
        if 'stacked' in file:
            continue
        ccd = CCDData.read(directory / file, unit="adu")
        ccd_list.append(ccd)

    stacked_bias = ccdproc.combine(ccd_list, method='median')
    stacked_bias.write(combined_bias_file, overwrite=True)


def make_main_bias_from_directory(directory: Path):
    df = get_list_of_files_in_directory(directory)
    # Filter for bias frames
    bias_frames = df[df['type'] == 'bias']

    # Group by CCD ID and binning
    grouped = bias_frames.groupby(['ccd_id', 'binning', 'read_mode'])

    # Produce stacked bias for each CCD and binning combination
    for (ccd_id, binning, read_mode), group in grouped:
        biases = group['filename'].tolist()
        if len(biases) > 0:
            stack_bias(directory, ccd_id, binning, read_mode, biases)


def make_main_flat_from_directory(directory: Path, flat_type: str):
    """

    Args:
        directory: Path to directory containing raw data
        flat_type: 'spec' or 'imag', whether we do flats of spectroscopic (long slit) or imaging frames

    Returns:

    """
    df = get_list_of_files_in_directory(directory)
    # Filter for flat frames
    flat_frames = df[df['type'] == 'flat']
    if flat_type == 'spec':
        flat_frames = flat_frames[flat_frames['slit_mask'].str.contains('ls')]
        # Group by CCD ID and slit mask and binning
        grouped = flat_frames.groupby(['ccd_id', 'binning', 'slit_mask'])
    else:
        flat_frames = flat_frames[flat_frames['slit_mask'].str.contains('imaging')]
        # Group by CCD ID and filter and binning
        grouped = flat_frames.groupby(['ccd_id', 'binning', 'filter'])

    # Produce stacked flats for each CCD, slit and binning combination
    for (ccd_id, binning, mask), group in grouped:
        flats = group['filename'].tolist()
        exposure_times = group['exptime'].tolist()
        if len(flats) > 0:
            if flat_type == 'spec':
                print(f'Combining {len(flats)} for ccd {ccd_id}, slit {mask} and binning {binning}')
                stack_flat(directory=directory, ccd_id=ccd_id, binning=binning,
                           files=flats, exposure_times=exposure_times, slit_mask=mask, photom_filter=None)
            else:
                print(f'Combining {len(flats)} for ccd {ccd_id}, filter {mask} and binning {binning}')
                stack_flat(directory=directory, ccd_id=ccd_id, binning=binning, files=flats,
                           exposure_times=exposure_times, slit_mask=None, photom_filter=mask)


def stack_flat(directory, ccd_id, binning, files, exposure_times, slit_mask=None, photom_filter=None, redo=False,
               method='median'):
    if slit_mask:
        combined_flat_file = directory / f'stacked_flat_{ccd_id}_{binning}_{slit_mask}.fits'
    elif photom_filter:
        # probably a sky flat.
        combined_flat_file = directory / f'stacked_flat_{ccd_id}_{binning}_{photom_filter}.fits'
    else:
        raise AssertionError("Flat must either have a slit mask or a photometric filter")

    if combined_flat_file.exists() and not redo:
        return
    ccd_list = []
    chosen_files = []
    for file, exptime in zip(files, exposure_times):
        if 'stacked' in file:
            continue
        ccd = CCDData.read(directory / file, unit="adu")
        ccd.data = ccd.data.astype(float)
        # remove the bias from the ccd, read from each file specifically
        read_mode = fits.getheader(directory / file).get('speed').lower()
        b = fits.getdata(directory / f'stacked_bias_{ccd_id}_{binning}_{read_mode}.fits')
        ccd.data -= b
        filtered_data = ccd.data[ccd.data > 0]
        median_value = np.nanpercentile(filtered_data, 90)
        print(f'{file} has 90th percentile value: {median_value:.01f}')
        if not (1000 < median_value < 50000):
            print('rejecting', file, f'(count value is {median_value:.0f})')
            continue
        ccd = ccd.divide(operand=median_value)
        ccd_list.append(ccd)
        chosen_files.append(file)
    if len(chosen_files) < 1:
        print("Not enough flat exposures for this configuration / CCD")
        return
    print(f"Combining files: {chosen_files}")

    stacked_flat = ccdproc.combine(ccd_list, method=method)
    # just making sure no zeros go through
    stacked_flat.data[stacked_flat.data < 1e-4] = np.nanmedian(stacked_flat.data)
    # and making the flat ~ 1
    flat_flat = stacked_flat.data.flatten()
    nonzero_flat = flat_flat[np.where(flat_flat > 0.)]
    norm = np.nanpercentile(nonzero_flat, 50)
    stacked_flat = stacked_flat.divide(norm)
    stacked_flat.write(combined_flat_file, overwrite=True)

