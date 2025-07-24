import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
import areeepo as ap
import niceplot
import h5py
import numpy as np
import argparse
from pathlib import Path
import cmasher as cmr
import astropy.units as au
import astropy.constants as ac
niceplot.set_bismuth_style(plt)

N = 128
def maxdens(f):
    dens = ap.load_gas(f, 'Density')
    units = ap.load_units(f)
    maxdens_iu = dens.max()
    iu2cgs = units['UnitDensity'] * au.g / au.cm**3
    iu2ncgs = (iu2cgs / ac.m_p).to(1./au.cm**3)
    print(iu2ncgs)
    maxdens_cgs = maxdens_iu * iu2cgs
    maxndens_cgs = (maxdens_cgs / ac.m_p).to(1./au.cm**3)

    return maxdens_iu, maxdens_cgs, maxndens_cgs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()

    for fname in args.files:
        with h5py.File(fname) as f:
            tmp = maxdens(f)
            for q in tmp:
                print(f"{q:.5e}  ", end="")
            print("")
    
    plt.show()
