import numpy as np

def load_two_col_csv(path: str):
    data = np.genfromtxt(path, delimiter=",", invalid_raise=False)

    # remove rows with NaNs (headers / corrupt lines)
    data = data[~np.isnan(data).any(axis=1)]

    wl = data[:, 0]
    y = data[:, 1]

    # sort by wavelength
    idx = np.argsort(wl)
    return wl[idx], y[idx]

