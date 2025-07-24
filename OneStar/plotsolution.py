import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
import areeepo as ap
import niceplot
import h5py
import numpy as np
import argparse
from pathlib import Path
import cmasher as cmr
style_file_path = Path(__file__).parent / 'plotstyle.mplstyle'
plt.style.use(style_file_path)
GFM_SOLAR_METALLICITY = 0.0127

def density_temperature_maps(fname, bins):
    with h5py.File(fname,'r') as f:
        pos = ap.load_gas(f,'Coordinates')
        dens = ap.load_gas(f,'Density')
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

def makeplot_maps(res, time_Myr, cmaps, labels, vmin=None, vmax=None):
    infotext = f"$t = {time_Myr:.2f}$ Myr"
    data = np.array(
        [
            [ np.log10(res[0]), np.log10(res[1]) ]
        ]
    )
    return niceplot.imshow_grid(1,2,data=data, s = 1.2, cmap=cmaps, label=labels,
                                infotext=infotext, infotext_fontsize=10, sharex=True, sharey=True,
                                vmin = vmin, vmax = vmax)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('bpath', nargs='+', type=str)
    parser.add_argument('--show', action='store_true')
    parser.add_argument('--vminmax', action='store_true')
    args = parser.parse_args()

    bins = [64,64,16]

    cmaps = [['cmr.ocean','magma']]
    vmin = [[-4.7 ,2.5 , 1e-9]]
    vmax = [[-4 ,4 , 4e-4]]
    labels = [[r'$\log_{10}\Sigma$ [g cm$^{-2}$]', r'$\log_{10}T$ [K]']]
    
    if not args.vminmax:
        vmin, vmax = None, None

    nmaps = 2
    nfiles = len(args.bpath)
    times = np.empty(nfiles)
    data = np.empty((nfiles, nmaps, bins[0], bins[1]))

    for i in range(nfiles):
        fpath = Path(args.bpath[i])
        res = density_temperature_maps(fpath, bins=bins)
        times[i] = res['time_Myr']
        data[i,0] = res['density']
        data[i,1] = res['temperature']
    
    bounds = np.zeros((nmaps,2))
    for k in range(nmaps):
        bounds[k,0] = data[:,k].min()
        bounds[k,1] = data[:,k].max()

    vmin = np.log10( bounds[:,0].reshape((1,nmaps)) )
    vmax = np.log10( bounds[:,1].reshape((1,nmaps)) )

    for i in range(nfiles):
        fig, ax, im, cbar = makeplot_maps(data[i], times[i], cmaps=cmaps, labels=labels, vmin=vmin, vmax=vmax)

        if args.show:
            plt.show()
        else:
            fname_out = fpath.parent / f"frame_{i:03n}.png"
            print(fname_out)
            fig.savefig(fname_out,dpi=200,bbox_inches='tight')

        plt.close(fig)