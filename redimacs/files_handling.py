import pandas as pd
from pathlib import Path
from astropy.io import fits


def get_list_of_files_in_directory(directory: Path):
    # List of FITS files in the directory
    fits_files = list(directory.glob('*.fits'))

    # DataFrame to store the metadata
    data = []

    # Iterate through each FITS file
    for fits_file in fits_files:
        # skip mosaics made by their iraf thingie
        if '_mos' in fits_file.name:
            continue
        with fits.open(fits_file) as hdul:
            header = hdul[0].header
            object_name = header.get('OBJECT', '')
            binning = header.get('BINNING', '1x1')
            if not object_name:
                continue
            object_name = object_name.lower()
            filename = fits_file.name
            ccd_id = 'c' + str(header.get('CHIP'))

            data.append({
                'filename': filename,
                'object': object_name,
                'ccd_id': ccd_id,
                'binning': binning,
                'type': header.get('EXPTYPE').lower(),
                'subraster': header.get('SUBRASTR'),
                'exptime': header.get('EXPTIME'),
                'slit_mask': header.get('SLITMASK', '').lower(),
                'filter': header.get('FILTER', '').lower()
            })

    df = pd.DataFrame(data)
    return df

