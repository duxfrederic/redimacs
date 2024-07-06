

def bin_spectrum(spectrum, bin_factor=3):
    num_rows, num_cols = spectrum.shape

    # Calculate the number of columns to keep after binning
    num_binned_cols = num_cols // bin_factor * bin_factor  # Ensure it's a multiple of bin_factor

    # Discard any extra columns
    spectrum_trimmed = spectrum[:, :num_binned_cols]

    # Reshape the array to perform binning
    spectrum_reshaped = spectrum_trimmed.reshape(num_rows, num_binned_cols // bin_factor, bin_factor)

    # Sum along the binning axis (axis 2)
    spectrum_binned = spectrum_reshaped.sum(axis=2)

    return spectrum_binned
