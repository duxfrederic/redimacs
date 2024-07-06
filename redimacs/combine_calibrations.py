import ccdproc
from astropy.nddata import CCDData
from astropy.io import fits
from pathlib import Path
import pandas as pd

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
    bias_frames = df[df['object'] == 'bias']

    # Group by CCD ID and binning
    grouped = bias_frames.groupby(['ccd_id', 'binning'])

    # Produce stacked bias for each CCD and binning combination
    for (ccd_id, binning), group in grouped:
        biases = group['filename'].tolist()
        if len(biases) > 0:
            print(f'Combining {len(biases)} for ccd {ccd_id} and binning {binning}')
            stack_bias(directory, ccd_id, binning, biases)
