################################################################################
# Copyright (c) 2024 Patrick Hirling (patrick.hirling@epfl.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# created by Patrick Hirling, last modified March 2025
################################################################################

import numpy as np
import h5py
import astropy.units as au
import argparse
from pathlib import Path
FloatType = np.float64
IntType = np.int32
from astropy.constants import G as G_Grav
import astropy.units as au
from scipy.spatial import KDTree

def imprint_turbulence_field(Pos, Vel, BoxSize, turbfile, target_sigma):
    with h5py.File(turbfile,'r') as fturb:
        turb_field_x = np.array(fturb['turb_field_x']).flatten()
        turb_field_y = np.array(fturb['turb_field_y']).flatten()
        turb_field_z = np.array(fturb['turb_field_z']).flatten()
        bins = np.array(fturb['N'], dtype=int)

    # Create cartesian grid
    X, dx = np.linspace(0,BoxSize,bins[0],retstep=True, endpoint=False)
    Y, dy = np.linspace(0,BoxSize,bins[1],retstep=True, endpoint=False)
    Z, dz = np.linspace(0,BoxSize,bins[2],retstep=True, endpoint=False)
    XX, YY, ZZ = np.meshgrid(X,Y,Z)
    gpoints = np.transpose(np.stack((XX.flatten(),YY.flatten(),ZZ.flatten())))

    # Create tree
    tree = KDTree(gpoints)

    # Find index of nearest neighbour to each grid point (voxel)
    dist, i = tree.query(Pos)
    
    Vel_out = np.copy(Vel)

    #print(Vel_out.shape, turb_field_x[i].shape)
    Vel_out[:, 0] = turb_field_x[i]
    Vel_out[:, 1] = turb_field_y[i]
    Vel_out[:, 2] = turb_field_z[i]

    comp_vel_disp = np.std(Vel_out, axis=0).sum()
    Vel_out *= target_sigma / comp_vel_disp
    return Vel_out

if __name__ == '__main__':
    import niceplot
    import areeepo as ap
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    turbfile = 'TurbGen_output.h5'
    ICfile = 'IC_SphericalCloudStatic_N64.hdf5'
    with h5py.File(ICfile, 'r') as f:
        Pos = np.array(f['PartType0/Coordinates'])
        Vel = np.array(f['PartType0/Velocities'])

    Vel2 = imprint_turbulence_field(Pos, Vel, 1.0,turbfile,22)
    VelNorm = np.linalg.norm(Vel2, axis=1)

    comp_vel_disp = np.std(Vel2, axis=0)
    print(comp_vel_disp)

    proj, bs, e2d = ap.project_with_NN(Pos, VelNorm, 2, 64)
    niceplot.imshow(proj)
    plt.show()