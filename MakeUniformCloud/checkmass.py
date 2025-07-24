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

def M_analytical(r, eps, fR):
    fact1 = (1 - eps)*4./3 * np.pi * fR**3 + eps
    invfact1 = 1.0 / fact1
    fact2 = 1 + eps*(1 + 3.0/(4*np.pi*fR**3))
    invfact2 = 1.0 / fact2
    rhoc = invfact1
    rhoa = eps * rhoc

    M1 = 4./3 * np.pi * r**3 * rhoc
    M2 = 4./3 * np.pi * (fR**3 * rhoc + r**3 * rhoa - fR**3 * rhoa)
    res = np.where(r <= fR, M1, M2)

    return plt.plot(r, res, color='black', lw=2, ls='--')

def make_mass_profile(f, rsp):
    pos = ap.load_gas(f, 'Coordinates')
    dens = ap.load_gas(f, 'Density')
    u = ap.load_gas(f, 'InternalEnergy')
    mass = ap.load_gas(f, 'Masses')
    bs = ap.get_boxsize(f)
    rad = np.linalg.norm(pos - bs/2., axis=1)
    nbins = len(rsp)
    res = np.empty(nbins)
    for i in range(nbins):
        idx = np.where(rad <= rsp[i])
        res[i] = mass[idx].sum()

    return plt.plot(rsp, res)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()

    nbins = 100
    rsp = np.linspace(0, 1, nbins)
    tmp2 = M_analytical(rsp, 1./100, 0.25)
    
    for fname in args.files:
        with h5py.File(fname) as f:
            tmp = make_mass_profile(f, rsp)


    plt.show()