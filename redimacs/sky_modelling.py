import numpy as np


def subtract_median(spectrum):
    num_rows, num_cols = spectrum.shape

    mask_start = num_rows // 4
    mask_end = num_rows - num_rows // 4

    mask = np.ones(num_rows, dtype=bool)
    mask[mask_start:mask_end] = False

    medians = np.nanmedian(spectrum[mask], axis=0)

    sky_subtracted = spectrum - medians

    return sky_subtracted
