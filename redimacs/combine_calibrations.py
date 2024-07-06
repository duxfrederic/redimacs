import ccdproc
from astropy.nddata import CCDData
from astropy.io import fits
from pathlib import Path
import pandas as pd
import numpy as np

from .files_handling import get_list_of_files_in_directory


def stack_bias(directory, ccd_id, binning, files, redo=False):
    combined_bias_file = directory / f'stacked_bias_{ccd_id}_{binning}.fits'
    if combined_bias_file.exists() and not redo:
        return
    ccd_list = []
    for file in files:
        if 'stacked' in file:
            continue
        ccd = CCDData.read(directory / file, unit="adu")
        ccd_list.append(ccd)

    stacked_bias = ccdproc.combine(ccd_list, method='average')
    stacked_bias.write(combined_bias_file, overwrite=True)


def make_main_bias_from_directory(directory: Path):
    df = get_list_of_files_in_directory(directory)
    # Filter for bias frames
    bias_frames = df[df['type'] == 'bias']

    # Group by CCD ID and binning
    grouped = bias_frames.groupby(['ccd_id', 'binning'])

    # Produce stacked bias for each CCD and binning combination
    for (ccd_id, binning), group in grouped:
        biases = group['filename'].tolist()
        if len(biases) > 0:
            print(f'Combining {len(biases)} for ccd {ccd_id} and binning {binning}')
            stack_bias(directory, ccd_id, binning, biases)


def make_main_flat_from_directory(directory: Path):
    df = get_list_of_files_in_directory(directory)
    # Filter for flat frames
    flat_frames = df[df['type'] == 'flat']
    flat_frames = flat_frames[flat_frames['slit_mask'].str.contains('ls')]

    # Group by CCD ID and slit mask and binning
    grouped = flat_frames.groupby(['ccd_id', 'binning', 'slit_mask'])

    # Produce stacked flats for each CCD, slit and binning combination
    for (ccd_id, binning, slit), group in grouped:
        flats = group['filename'].tolist()
        exposure_times = group['exptime'].tolist()
        if len(flats) > 0:
            print(f'Combining {len(flats)} for ccd {ccd_id}, slit {slit} and binning {binning}')
            stack_flat(directory, ccd_id, binning, slit, flats, exposure_times)


def stack_flat(directory, ccd_id, binning, slit_mask, files, exposure_times, redo=False):
    combined_flat_file = directory / f'stacked_flat_{ccd_id}_{binning}_{slit_mask}.fits'
    if combined_flat_file.exists() and not redo:
        return
    ccd_list = []
    for file, exptime in zip(files, exposure_times):
        if 'stacked' in file:
            continue
        ccd = CCDData.read(directory / file, unit="adu")
        large_value = np.nanpercentile(ccd.data, 95)
        if not (5000 < large_value < 40000):
            print('rejecting', file, f'(count value is {large_value:.0f}')
        ccd.divide(operand=exptime)
        ccd_list.append(ccd)

    stacked_flat = ccdproc.combine(ccd_list, method='median')
    stacked_flat.divide(np.nanpercentile(stacked_flat.data, 50))
    stacked_flat.write(combined_flat_file, overwrite=True)

