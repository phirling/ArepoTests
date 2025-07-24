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

def density_temperature_maps(fname, bins):
    with h5py.File(fname,'r') as f:
        pos = ap.load_gas(f,'Coordinates')
        dens = ap.load_gas(f,'Density')
        #temp = ap.load_gas(f,'InternalEnergy')
        temp = ap.compute_temperature(f)
        units = ap.load_units(f)
        is_cosmo = ap.is_cosmological(f)
        time = ap.get_time(f)
        if is_cosmo:
            zred = 1./time - 1
            time_Myr = ((Planck15.age(zred) - Planck15.age(9.0)).to('Myr')).value
        else:
            time_Myr = time * units['UnitTime'] / ap.Myr_cgs
    fields = [dens, temp*dens]
    proj, bs, e2d = ap.project_with_NN(pos, fields, 2, bins)

    densfact = units['UnitMass'] / units['UnitLength']**2
    if is_cosmo:
        densfact *= (1./time)**2
    res = {'binsizes' : bs, 'extent' : e2d, 'time_Myr' : time_Myr, 'units' : units}
    res['density'] = proj[0] * densfact
    res['temperature'] = proj[1] / proj[0]

    return res

def get_stars(fname):
    with h5py.File(fname,'r') as f:
        if 'PartType4' in f.keys():
            spos = ap.load_dataset_from_parttype(f, 4, 'Coordinates', remove_h_factors=1)
            return spos
        else:
            return None

def makeplot_maps(vals, time_Myr, cmaps, labels, vmin=None, vmax=None, extent = None):
    infotext = f"$t = {time_Myr:.2f}$ Myr"
    data = np.array(
        [
            [ np.log10(vals[0]), np.log10(vals[1]) ]
        ]
    )
    return niceplot.imshow_grid(1,2,data=data, s = 1.6, cmap=cmaps, label=labels, extent = extent,
                                infotext=infotext, infotext_fontsize=10, sharex=True, sharey=True,
                                vmin = vmin, vmax = vmax)

def get_index_of_PID(fname, PID):
    with h5py.File(fname,'r') as f:
        if 'PartType4' in f.keys():
            fPID = ap.load_dataset_from_parttype(f, 4, 'ParticleIDs', remove_h_factors=0)
            return int(np.where(fPID == PID)[0][0])
        else:
            return None
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('bpath', nargs='+', type=str)
    parser.add_argument('-N', type=int, default=64)
    parser.add_argument('-Nz', type=int, default=16)
    parser.add_argument('--show', action='store_true')
    parser.add_argument('--vminmax', action='store_true')
    args = parser.parse_args()

    bins = [int(args.N),int(args.N),int(args.Nz)]

    cmaps = [['cmr.ocean','magma']]
    labels = [[r'$\log_{10}\Sigma$ [g cm$^{-2}$]', r'$\log_{10}T$ [K]']]
    
    vmin, vmax = None, None
    PID_select = 10441 # Found by hand

    nmaps = 2
    nfiles = len(args.bpath)
    times = np.empty(nfiles)
    data = np.empty((nfiles, nmaps, bins[0], bins[1]))
    star_positions = []
    star_indices = []
    extent = None
    print("Making projections...")
    for i in tqdm(range(nfiles)):
        fpath = Path(args.bpath[i])
        res = density_temperature_maps(fpath, bins=bins)
        times[i] = res['time_Myr']
        data[i,0] = res['density']
        data[i,1] = res['temperature']
        star_positions.append(get_stars(fpath))
        star_indices.append(get_index_of_PID(fpath, PID_select))
        if i == 0: extent = res['extent']

    bounds = np.zeros((nmaps,2))
    for k in range(nmaps):
        bounds[k,0] = data[:,k].min()
        bounds[k,1] = data[:,k].max()

    vmin = np.log10( bounds[:,0].reshape((1,nmaps)) )
    vmax = np.log10( bounds[:,1].reshape((1,nmaps)) )

    print("Rendering images...")
    for i in range(nfiles):
        fig, ax, im, cbar = makeplot_maps(data[i], times[i], cmaps=cmaps, labels=labels, vmin=vmin, vmax=vmax, extent=extent)
        if star_positions[i] is not None:
            ax[0,0].scatter(star_positions[i][:,0], star_positions[i][:,1], s=4, c='paleturquoise', marker='*')
        if star_indices[i] is not None:
            j = star_indices[i]
            ax[0,0].scatter(star_positions[i][j,0], star_positions[i][j,1], s=150, c='orange', ec='black', marker='*')
        for axx in ax.flatten():
            axx.set_xlim(extent[0:2])
            axx.set_ylim(extent[2:4])

        if args.show:
            plt.show()
        else:
            fname_out = fpath.parent / f"frame_{i:03n}.png"
            print(fname_out)
            fig.savefig(fname_out,dpi=200,bbox_inches='tight')

        plt.close(fig)
    print("done.")