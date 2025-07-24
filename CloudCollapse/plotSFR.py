import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
import areeepo as ap
import niceplot
import h5py
import numpy as np
import argparse
from pathlib import Path
import cmasher as cmr
from tqdm import tqdm
niceplot.set_bismuth_style(plt)
GFM_SOLAR_METALLICITY = 0.0127

def get_stellar_mass(fname):
    with h5py.File(fname,'r') as f:
        if 'PartType4' in f.keys():
            smass = ap.load_dataset_from_parttype(f, 4, 'Masses', remove_h_factors=1)

            return smass.sum()
        else:
            return 0.0

def get_gas_mass(fname):
    with h5py.File(fname,'r') as f:
        gmass = ap.load_dataset_from_parttype(f, 0, 'Masses', remove_h_factors=1)
    return gmass.sum()

def get_time(fname):
    with h5py.File(fname,'r') as f:
        is_cosmo = ap.is_cosmological(f)
        units = ap.load_units(f)
        time = ap.get_time(f)
        if is_cosmo:
            zred = 1./time - 1
            time_Myr = ((Planck15.age(zred) - Planck15.age(9.0)).to('Myr')).value
        else:
            time_Myr = time * units['UnitTime'] / ap.Myr_cgs
        return time_Myr
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('bpath', nargs='+', type=str)
    parser.add_argument('--show', action='store_true')
    args = parser.parse_args()

    nfiles = len(args.bpath)
    times = np.empty(nfiles)
    stellar_masses = np.empty(nfiles)
    gas_masses = np.empty(nfiles)

    for i in range(nfiles):

        fpath = Path(args.bpath[i])
        stellar_masses[i] = get_stellar_mass(fpath)
        gas_masses[i] = get_gas_mass(fpath)
        times[i] = get_time(fpath)

    sfr = np.diff(stellar_masses) / np.diff(times) / 1e6 # Msun / yr

    fig, ax = plt.subplots(2,2, figsize=(10,8))
    ax[0,0].plot(times, stellar_masses/1e6)
    ax[0,0].set_xlabel('Time [Myr]')
    ax[0,0].set_ylabel(r'Stellar Mass [$10^6$ M$_\odot$]')
    ax[0,1].plot(times[:-1], sfr)
    ax[0,1].set_xlabel('Time [Myr]')
    ax[0,1].set_ylabel(r'SFR [M$_\odot$ yr$^{-1}$]')

    ax[1,0].plot(times, stellar_masses/1e6, label = 'Stars')
    ax[1,0].plot(times, gas_masses/1e6, label = 'Gas')
    ax[1,0].set_xlabel('Time [Myr]')
    ax[1,0].set_ylabel(r'Mass [$10^6$ M$_\odot$]')
    ax[1,0].legend()

    ax[1,1].plot(times, stellar_masses / gas_masses[0])
    ax[1,1].set_xlabel('Time [Myr]')
    ax[1,1].set_ylabel(r'Star Formation Efficiency')

    if args.show:
        plt.show()
    else:
        fname_out = "SFR.png"
        fig.savefig(fname_out,dpi=200,bbox_inches='tight')

        plt.close(fig)
    print("done.")



"""
    fig, ax = plt.subplots(1,2, figsize=(12,4))
    ax[0].plot(times, stellar_masses/1e6)
    ax[0].set_xlabel('Time [Myr]')
    ax[0].set_ylabel(r'Stellar Mass [$10^6$ M$_\odot$]')
    ax[1].plot(times[:-1], sfr)
    ax[1].set_xlabel('Time [Myr]')
    ax[1].set_ylabel(r'SFR [M$_\odot$ yr$^{-1}$]')

    fig2, ax2 = plt.subplots(1,2, figsize=(12,4))
    ax2[0].plot(times, stellar_masses/1e6, label = 'Stars')
    ax2[0].plot(times, gas_masses/1e6, label = 'Gas')
    ax2[0].set_xlabel('Time [Myr]')
    ax2[0].set_ylabel(r'Mass [$10^6$ M$_\odot$]')
    ax2[0].legend()

    ax2[1].plot(times, stellar_masses / gas_masses[0])
    ax2[1].set_xlabel('Time [Myr]')
    ax2[1].set_ylabel(r'Star Formation Efficiency')
    """