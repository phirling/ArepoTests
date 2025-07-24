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
from plotmetals import GFM_N_CHEM_ELEMENTS, compute_global_metal_masses, get_time, element_labels, get_star_mass

def sname(dir, i):
    return dir + f"snap_{i:03n}.hdf5"

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=str, default="")
    parser.add_argument('--suffix', type=str, default="")
    parser.add_argument('-mode', type=int, default=0)
    parser.add_argument('-o', type=str, default=None)
    args = parser.parse_args()

    mode = int(args.mode)
    prefix = str(args.prefix)
    suffix = str(args.suffix)

    nfiles = 504
    dir_full = prefix + 'full/' + suffix
    dir_SNII = prefix + 'onlySNII/' + suffix
    dir_SNIa = prefix + 'onlySNIa/' + suffix
    dir_AGB  = prefix + 'onlyAGB/' + suffix

    data_elements = np.empty((4, nfiles, GFM_N_CHEM_ELEMENTS))
    times = np.empty(nfiles)

    initial_star_mass = get_star_mass(sname(dir_full, 0))

    for i in range(nfiles):
        times[i] = get_time(sname(dir_full, i))

        data_elements[0,i] = compute_global_metal_masses(sname(dir_full, i))
        data_elements[1,i] = compute_global_metal_masses(sname(dir_SNII, i))
        data_elements[2,i] = compute_global_metal_masses(sname(dir_SNIa, i))
        data_elements[3,i] = compute_global_metal_masses(sname(dir_AGB, i))

    s0, s1, s2, s3 = np.sum(data_elements[0, :,2:], axis=1), np.sum(data_elements[1, :,2:], axis=1), np.sum(data_elements[2, :,2:], axis=1), np.sum(data_elements[3, :,2:], axis=1)

    cmap = plt.get_cmap('rainbow')
    s = 1.2
    fig, ax = plt.subplots(figsize=(s*6,s*6))

    if mode == 0 or mode == 1:
        if mode == 0:
            itr = range(GFM_N_CHEM_ELEMENTS-2)
        else:
            itr = [0, 2, 6]
        for k in itr:
            c = cmap(k / (GFM_N_CHEM_ELEMENTS-3))
            ax.semilogy(times, data_elements[0, :, k+2], color = c, ls='-', label=element_labels[k+2])
            ax.semilogy(times, data_elements[1, :, k+2], color = c, ls='--')
            ax.semilogy(times, data_elements[2, :, k+2], color = c, ls='-.')
            ax.semilogy(times, data_elements[3, :, k+2], color = c, ls=':')

        lw=2.5

        ax.semilogy(times, s0, color='black', ls='-',  lw=lw, label = 'Full')
        ax.semilogy(times, s1, color='black', ls='--', lw=lw, label = 'SNII')
        ax.semilogy(times, s2, color='black', ls='-.',  lw=lw, label = 'SNIa')
        ax.semilogy(times, s3, color='black', ls=':', lw=lw, label = 'AGB')
        
        if mode == 0:
            ax.legend(loc='lower center', bbox_to_anchor = (0.5,1), ncols=6)
        else:
            ax.legend(loc='lower center', bbox_to_anchor = (0.5,1), ncols=4)

    elif mode == 2:
        lw=2
        ax.semilogy(times, s0, color='black', ls='-',  lw=lw, label = 'Full')
        ax.semilogy(times, s1, color='black', ls='--', lw=lw, label = 'SNII')
        ax.semilogy(times, s2, color='black', ls='-.',  lw=lw, label = 'SNIa')
        ax.semilogy(times, s3, color='black', ls=':', lw=lw, label = 'AGB')

        ax.semilogy(times, s1 + s2 + s3, color='red', ls='--', lw=lw*1.5, label = 'Sum of Individuals')
        ax.legend(loc='lower center', bbox_to_anchor = (0.5,1), ncols=3)
    
    elif mode == 3:
        itr = range(GFM_N_CHEM_ELEMENTS)
        for k in itr:
            c = cmap(k / (GFM_N_CHEM_ELEMENTS-1))
            ax.semilogy(times, data_elements[0, :, k] / data_elements[0, 0, k], color = c, ls='-', label=element_labels[k])
            ax.semilogy(times, data_elements[1, :, k] / data_elements[1, 0, k], color = c, ls='--')
            ax.semilogy(times, data_elements[2, :, k] / data_elements[2, 0, k], color = c, ls='-.')
            ax.semilogy(times, data_elements[3, :, k] / data_elements[3, 0, k], color = c, ls=':')

        lw=2.5

        ax.semilogy(times, s0, color='black', ls='-',  lw=lw, label = 'Full')
        ax.semilogy(times, s1, color='black', ls='--', lw=lw, label = 'SNII')
        ax.semilogy(times, s2, color='black', ls='-.',  lw=lw, label = 'SNIa')
        ax.semilogy(times, s3, color='black', ls=':', lw=lw, label = 'AGB')
        ax.legend(loc='lower center', bbox_to_anchor = (0.5,1), ncols=6)

    elif mode == 4:
        n = 2
        itr = range(n)
        for k in itr:
            c = cmap(k / (n-1))
            ax.semilogy(times, (data_elements[0, :, k] - data_elements[0, 0, k]) / data_elements[0, 0, k] / 100, color = c, ls='-', label=element_labels[k])
            ax.semilogy(times, (data_elements[1, :, k] - data_elements[1, 0, k]) / data_elements[1, 0, k] / 100, color = c, ls='--')
            ax.semilogy(times, (data_elements[2, :, k] - data_elements[2, 0, k]) / data_elements[2, 0, k] / 100, color = c, ls='-.')
            ax.semilogy(times, (data_elements[3, :, k] - data_elements[3, 0, k]) / data_elements[3, 0, k] / 100, color = c, ls=':')

        lw=2.5
        ax.legend(loc='lower center', bbox_to_anchor = (0.5,1), ncols=6)
        ax.set_ylabel(r"Rel. Delta Mass [%]")

    ax.set_xlabel("Time [Myr]")
    if mode not in [4]:
        ax.set_ylabel(r"Mass [$\mathrm{M_\odot}$]")
    fig.tight_layout()
    if args.o is None:
        plt.show()
    else:
        fig.savefig(str(args.o), dpi=300, bbox_inches='tight')