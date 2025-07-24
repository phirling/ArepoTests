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
    parser.add_argument('-o', type=str, default=None)
    parser.add_argument('-i', type=int, default=100)
    args = parser.parse_args()
    
    prefix = str(args.prefix)
    suffix = str(args.suffix)

    i_initial_snap = 0
    i_final_snap = int(args.i)
    dir_full = prefix + 'full/' + suffix
    dir_SNII = prefix + 'onlySNII/' + suffix
    dir_SNIa = prefix + 'onlySNIa/' + suffix
    dir_AGB  = prefix + 'onlyAGB/' + suffix

    data_elements_init = np.empty((4, GFM_N_CHEM_ELEMENTS))
    data_elements_final = np.empty((4, GFM_N_CHEM_ELEMENTS))

    initial_star_mass = get_star_mass(sname(dir_full, i_initial_snap))

    time_final = get_time(sname(dir_full, i_final_snap))

    data_elements_init[0] = compute_global_metal_masses(sname(dir_full, i_initial_snap))
    data_elements_init[1] = compute_global_metal_masses(sname(dir_SNII, i_initial_snap))
    data_elements_init[2] = compute_global_metal_masses(sname(dir_SNIa, i_initial_snap))
    data_elements_init[3] = compute_global_metal_masses(sname(dir_AGB, i_initial_snap))

    data_elements_final[0] = compute_global_metal_masses(sname(dir_full, i_final_snap))
    data_elements_final[1] = compute_global_metal_masses(sname(dir_SNII, i_final_snap))
    data_elements_final[2] = compute_global_metal_masses(sname(dir_SNIa, i_final_snap))
    data_elements_final[3] = compute_global_metal_masses(sname(dir_AGB, i_final_snap))

    #sums = np.sum(data_elements_final, axis=1)
    mass_frac = (data_elements_final - data_elements_init) / initial_star_mass
    #mass_frac_tot = sums / initial_star_mass
    

    cmap = plt.get_cmap('rainbow')
    s = 1.2

    fig, ax = plt.subplots(figsize=(s*6,s*4))
    ax.set_yscale('log')
    lbls = element_labels[:-1]
    # 0.7 1.7
    x1 = 0.7 + np.arange(len(lbls))
    x2 = 0.9 + np.arange(len(lbls))
    x3 = 1.1 + np.arange(len(lbls))
    x4 = 1.3 + np.arange(len(lbls))
    xticks = 1.0 + np.arange(len(lbls))
    width = 0.2

    mass_frac[2,0] = 0
    mass_frac[2,1] = 0
    mass_frac[2,3] = 0

    ax.bar(x1, mass_frac[0,:-1], width=width, label = 'Total')
    ax.bar(x2, mass_frac[1,:-1], width=width, label = 'SNII')
    ax.bar(x3, mass_frac[3,:-1], width=width, label = 'AGB')
    ax.bar(x4, mass_frac[2,:-1], width=width, label = 'SNIa')

    ax.set_xticks(xticks, lbls)
    ax.set_ylabel("Returned stellar mass fraction")
    ax.legend()
    ax.set_title(f"Time = {time_final/1000:.2f} Gyr")
    ax.set_ylim(2e-6,0.4)
    #ax.bar(element_labels[:-1], mass_frac[1,:-1])
    if args.o is None:
        plt.show()
    else:
        fig.savefig(str(args.o), dpi=300, bbox_inches='tight')
    
    #ax.set_xlabel("Time [Myr]")
    #ax.set_ylabel(r"Mass [$\mathrm{M_\odot}$]")
    #fig.tight_layout()
    #plt.show()