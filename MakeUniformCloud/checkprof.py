import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
import areeepo as ap
import niceplot
import h5py
import numpy as np
import argparse
from pathlib import Path
import cmasher as cmr
niceplot.set_bismuth_style(plt)

N = 128
def make_profile(f):
    pos = ap.load_gas(f, 'Coordinates')
    dens = ap.load_gas(f, 'Density')
    u = ap.load_gas(f, 'InternalEnergy')
    bs = ap.get_boxsize(f)
    rad = np.linalg.norm(pos - bs/2., axis=1)

    return plt.scatter(rad, dens, s=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()

    fname_ic = str(args.files[0])
    fname_relax = str(args.files[1])

    for fname in args.files:
        with h5py.File(fname) as f:
            tmp = make_profile(f)
    
    plt.show()