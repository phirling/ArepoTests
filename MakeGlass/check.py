import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
import areeepo as ap
import niceplot
import h5py
import numpy as np
import argparse
from pathlib import Path
import cmasher as cmr

N = 64
def make_projections(f_ic, f_relax):
    pos1 = ap.load_gas(f_ic, 'Coordinates')
    dens1 = ap.load_gas(f_ic, 'Density')
    u1 = ap.load_gas(f_ic, 'InternalEnergy')
    quant1 = [dens1, u1*dens1]

    bs = ap.get_boxsize(f_ic)
    extent = np.array([
        [0,bs],
        [0,bs],
        [0,bs]
    ])
    res1, _, _ = ap.project_with_NN(pos1, quant1, 2, bins=N, extent=extent)

    pos2 = ap.load_gas(f_relax, 'Coordinates')
    dens2 = ap.load_gas(f_relax, 'Density')
    u2 = ap.load_gas(f_relax, 'InternalEnergy')
    quant2 = [dens2, u2*dens2]

    res2, _, _ = ap.project_with_NN(pos2, quant2, 2, bins=N, extent=extent)

    data = np.array([
        [res1[0], res1[1] / res1[0]],
        [res2[0], res2[1] / res2[0]],
    ])

    m1, m2, m3, m4 = np.mean(res1[0]), np.mean(res1[1]/res1[0]), np.mean(res2[0]), np.mean(res2[1]/res2[0])
    data = np.array([
        [res1[0]/m1, res1[1] / res1[0]/m2],
        [res2[0]/m3, res2[1] / res2[0]/m4],
    ])
    return niceplot.imshow_grid(2,2,data, cmap='RdBu')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs=2)
    args = parser.parse_args()

    fname_ic = str(args.files[0])
    fname_relax = str(args.files[1])

    with h5py.File(fname_ic) as f1, h5py.File(fname_relax) as f2:
        tmp = make_projections(f1,f2)
    
    plt.show()