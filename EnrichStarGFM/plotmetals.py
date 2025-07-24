import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
import areeepo as ap
import niceplot
import h5py
import numpy as np
import argparse
from pathlib import Path
import cmasher as cmr
from warnings import warn
style_file_path = Path(__file__).parent / 'plotstyle.mplstyle'
plt.style.use(style_file_path)
GFM_N_CHEM_ELEMENTS = 10

# H, He, C, N, O, Ne, Mg, Si and Fe
element_labels = [
    'H',
    'He',
    'C',
    'N',
    'O',
    'Ne',
    'Mg',
    'Si',
    'Fe',
    'Others'
]

def compute_total_mass(fname):
    with h5py.File(fname,'r') as f:
        mass = ap.load_gas(f, 'Masses', remove_h_factors=1)
    total_mass = mass.sum()
    return total_mass

def get_star_mass(fname):
    with h5py.File(fname,'r') as f:
        smass = ap.load_dataset_from_parttype(f, 4, 'Masses', remove_h_factors=1)
    if len(smass) != 1:
        warn("More than 1 star found in file !")
    return smass[0]

def compute_global_metallicity(fname):
    with h5py.File(fname,'r') as f:
        Z = ap.load_gas(f, 'GFM_Metallicity', remove_h_factors=0)
        mass = ap.load_gas(f, 'Masses')

    total_mass = mass.sum()
    total_metal_mass = np.sum(Z * mass)

    return total_metal_mass / total_mass

def compute_global_metal_mass_from_metallicity(fname):
    with h5py.File(fname,'r') as f:
        Z = ap.load_gas(f, 'GFM_Metallicity', remove_h_factors=0)
        mass = ap.load_gas(f, 'Masses', remove_h_factors=1)

    total_metal_mass = np.sum(Z * mass)

    return total_metal_mass

def compute_global_metal_fractions(fname):
    with h5py.File(fname,'r') as f:
        Metals = ap.load_gas(f, 'GFM_Metals', remove_h_factors=0)
        mass = ap.load_gas(f, 'Masses', remove_h_factors=1)

    total_mass = mass.sum()
    element_masses = np.einsum('i...,i', Metals, mass)

    return element_masses / total_mass

def compute_global_metal_masses(fname):
    with h5py.File(fname,'r') as f:
        Metals = ap.load_gas(f, 'GFM_Metals', remove_h_factors=0)
        mass = ap.load_gas(f, 'Masses', remove_h_factors=1)

    element_masses = np.einsum('i...,i', Metals, mass)

    return element_masses

def get_time(fname):
    with h5py.File(fname,'r') as f:
        time = ap.get_time(f)
        is_cosmo = ap.is_cosmological(f)
        if is_cosmo:
            zred = 1./time - 1
            time_Myr = ((Planck15.age(zred) - Planck15.age(9.0)).to('Myr')).value
        else:
            units = ap.load_units(f)
            time_Myr = time * units['UnitTime'] / ap.Myr_cgs
    
    return time_Myr


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('bpath', nargs='+', type=str)
    parser.add_argument('--show', action='store_true')
    parser.add_argument('-o', type=str, default=None)
    args = parser.parse_args()


    nfiles = len(args.bpath)
    times = np.empty(nfiles)
    data = np.empty((nfiles, 1 + GFM_N_CHEM_ELEMENTS))
    data2 = np.empty((nfiles, 1))

    for i in range(nfiles):
        fpath = Path(args.bpath[i])
        data[i,0] = compute_global_metal_mass_from_metallicity(fpath)
        times[i] = get_time(fpath)
        metals = compute_global_metal_masses(fpath)
        data[i, 1:] = metals
        data2[i, 0] = get_star_mass(fpath)

    fig, ax = plt.subplots()

    cmap = plt.get_cmap('rainbow')
    for k in range(GFM_N_CHEM_ELEMENTS-2):
        c = cmap(k / (GFM_N_CHEM_ELEMENTS-3))
        ax.semilogy(times, data[:, k+3], color = c, ls='-', label=element_labels[k+2])
    ax.semilogy(times, np.sum(data[:,3:], axis=1), color='red', ls='-', lw=2)
    ax.semilogy(times, data[:,0], color='black', ls='--', lw=2)

    #ax.set_ylim(1e-7, None)
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel(r"Mass [$\mathrm{M_\odot}$]")

    ax.legend(ncols=2)

    fig2, ax2 = plt.subplots()
    ax2.plot(times, data2[:,0] / data2[0,0])
    ax2.set_xlabel("Time [Myr]")
    ax2.set_ylabel(r"$M_\star / M_\star^0$")
    if args.show:
        plt.show()
    if args.o is not None:
        fig.savefig(str(args.o), dpi=300, bbox_inches='tight')

    plt.close(fig)