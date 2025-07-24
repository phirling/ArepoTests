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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', type=str, default='static', choices=['static','gaussian','turb'])
    args = parser.parse_args()

    mode = str(args.mode)

    # Main physical parameters
    CLOUD_RADIUS_PC = 10.0
    CLOUD_MASS_SOLAR = 1e6
    EPS = 1. / 100
    FR = 1. / 4
    GFM_N_CHEM_ELEMENTS = 10
    GFM_SOLAR_METALLICITY = 0.0127
    BASE_IC_PATH = Path('/Users/patrickhirling/Galspec/misc/ICs/SphericalCloud_N64.hdf5')
    AMBIENT_TEMPERATURE = 1e4 # in Kelvin. Cloud temperature is a factor EPS lower
    MU = 1.26 # (We assume an ionized primordial gas here, approx)
    GAMMA = 5.0 / 3.0
    ALPHAVIR = 10
    IC_OUTPUT_FILENAME = 'IC_SphericalCloudStatic_N64.hdf5'

    # -------------------- Constants
    m_p_cgs = 1.67262158e-24 # g
    k_b_cgs = 1.380649e-16
    G_cgs = 6.6743e-08 # cm3 / (g s2)
    Msun_cgs = 1.9891e33
    HubbleParam = 0.6774

    # -------------------- Units
    UnitLength_in_cm            = 3.08568e21    # 1 kpc
    UnitMass_in_g               = Msun_cgs     # 1 Solar Mass
    UnitVelocity_in_cm_per_s    = 100000        # km / s
    UnitTime_in_s               = UnitLength_in_cm / UnitVelocity_in_cm_per_s
    UnitDensity_in_g_per_cm3 = UnitMass_in_g / UnitLength_in_cm**3
    UnitEnergy = UnitMass_in_g * UnitLength_in_cm**2 / UnitTime_in_s**2
    uts = UnitTime_in_s * au.s

    # -------------------- Set Main Parameters
    BoxSize = FloatType(CLOUD_RADIUS_PC / FR / 1000.0)  # BoxSize is in internal units (kpc)
    # For cosmological run
    TimeBegin = 0.1

    # -------------------- Setup Gas
    # Load Glass Cube ICs
    with h5py.File(BASE_IC_PATH) as f:
        Pos = np.array(f['PartType0/Coordinates'])
        Mass = np.array(f['PartType0/Masses'])
        Npart_gas = len(Pos)
        Rad = np.linalg.norm(Pos - 0.5, axis=1)
        idx_cloud = np.where(Rad <= FR)

    # Rescale the glass ICs to desired values
    fact1 = (1 - EPS)*4./3 * np.pi * FR**3 + EPS
    invfact1 = 1.0 / fact1
    fact2 = 1 + EPS*(1 + 3.0/(4*np.pi*FR**3))
    invfact2 = 1.0 / fact2
    
    mass_fact = CLOUD_MASS_SOLAR * fact2
    target_total_mass = mass_fact
    Pos *= BoxSize
    Mass *= mass_fact

    # Set internal energy
    cell_u = ( k_b_cgs*AMBIENT_TEMPERATURE / ((GAMMA-1) * MU * m_p_cgs) ) / UnitEnergy * UnitMass_in_g
    Utherm = np.full(Npart_gas, cell_u, dtype=FloatType)
    Utherm[idx_cloud] *= EPS

    # Velocities
    # Virial parameter
    veldisp_raw = np.sqrt(3 * G_Grav * ALPHAVIR * CLOUD_MASS_SOLAR*au.Msun / (5*CLOUD_RADIUS_PC*au.pc))
    veldisp = veldisp_raw.to(au.km/au.s)
    veldisp_per_dim = veldisp.value / 3

    Vel = np.zeros((Npart_gas, 3), dtype=FloatType)

    if mode == 'gaussian':
        Npart_cloud = len(idx_cloud)
        Vel_rand = np.random.normal(loc=0.0, scale=veldisp_per_dim, size=(Npart_gas, 3))
        Vel[idx_cloud] = Vel_rand[idx_cloud]

    elif mode == 'turb':
        from imprint_turbulence import imprint_turbulence_field
        Vel_turb = imprint_turbulence_field(Pos, Vel, BoxSize, 'TurbGen_output.h5', veldisp.value)
        Vel[idx_cloud] = Vel_turb[idx_cloud]

    compstd = np.std(Vel[idx_cloud], axis=0)
    print(f"computed velocity divergence (IU):", compstd)

    # -------------------- Additional Fields for SGChem
    # Abundances (order: CarbAbund, OxyAbund, MAbund, ZAtom)
    # We add a "metallicity floor" to avoid negative abundances due to advection
    Abund = np.zeros((Npart_gas, 4), dtype=FloatType)
    metallicity_floor_massfrac = 1e-10
    # Factors to convert from mass fraction to number fraction (abundance)
    factor_Carbon = 1. / (12 * 0.76)
    factor_Oxygen= 1. / (16 * 0.76)
    Abund[:,0] = metallicity_floor_massfrac * factor_Carbon
    Abund[:,1] = metallicity_floor_massfrac * factor_Oxygen
    Abund[:,2] = 1e-09 # Don't care about M
    Abund[:,3] = 8 * metallicity_floor_massfrac / GFM_SOLAR_METALLICITY # Sum of all metals in GFM

    # Dust to gas ratio (set to 0.1 by default ?)
    DustRatio = np.full(Npart_gas, 0.1, dtype=FloatType)

    # -------------------- Apply cosmological factors
    a = TimeBegin
    h = HubbleParam

    # ============================= Create Output File ============================= #
    print("Writing file: ", IC_OUTPUT_FILENAME, "...")
    zred_init = 1. / TimeBegin - 1
    Redshift = FloatType(zred_init)

    with h5py.File(IC_OUTPUT_FILENAME, 'w') as f:
        header = f.create_group('Header')
        gp_gas = f.create_group('PartType0')
        gp_star = f.create_group('PartType4')

        # Header
        NumPart = np.array([Npart_gas, 0, 0, 0, 0, 0], dtype=IntType)
        header.attrs.create('NumPart_ThisFile', NumPart)
        header.attrs.create('NumPart_Total', NumPart)
        header.attrs.create('NumPart_Total_HighWord', np.zeros(6, dtype=IntType))
        header.attrs.create('MassTable', np.zeros(6, dtype=IntType))
        header.attrs.create('Time', TimeBegin)
        header.attrs.create('Redshift', Redshift)
        header.attrs.create('BoxSize', BoxSize / a * h)
        header.attrs.create('NumFilesPerSnapshot', 1)
        header.attrs.create('Omega0', 0.0)
        header.attrs.create('OmegaB', 0.0)
        header.attrs.create('OmegaLambda', 0.0)
        header.attrs.create('HubbleParam', 1.0)
        header.attrs.create('Flag_Sfr', 0)
        header.attrs.create('Flag_Cooling', 0)
        header.attrs.create('Flag_StellarAge', 0)
        header.attrs.create('Flag_Metals', 0)
        header.attrs.create('Flag_Feedback', 0)

        if Pos.dtype == np.float64:
            header.attrs.create('Flag_DoublePrecision', 1)
        else:
            header.attrs.create('Flag_DoublePrecision', 0)

        # PartType0 Dataset (Gas)
        gp_gas.create_dataset('ParticleIDs', data = np.arange(1, Npart_gas + 1))
        gp_gas.create_dataset('Coordinates', data = Pos / a * h)
        gp_gas.create_dataset('Masses', data = Mass * h)
        gp_gas.create_dataset('Velocities', data = Vel / a**(0.5))
        gp_gas.create_dataset('InternalEnergy', data = Utherm)
        gp_gas.create_dataset('ElementAbundances', data = Abund)
        gp_gas.create_dataset('DusttoGasRatio', data = DustRatio)

    rho_mean_cloud = CLOUD_MASS_SOLAR*au.Msun / (4./3 * np.pi * (CLOUD_RADIUS_PC*au.pc)**3)
    freefall_time = np.sqrt(3*np.pi / (32*G_Grav*rho_mean_cloud))
    print(f"Estimated TargetGasMass: {target_total_mass/Npart_gas * h:.3g}. Mean effective cell mass: {(Mass*h).mean():.3g}")
    print(f"BoxSize: {BoxSize} (physical kpc) = {BoxSize / a * h} (comoving kpc, times h)")
    print(f"Velocity dispersion for alpha_vir={ALPHAVIR:.2f} : {veldisp:.3f}")
    print(f"Unit Time = {uts.to('Gyr'):.3g}")
    print(f"Free-fall time = {freefall_time.to('Myr'):.4g}")
    print("done.")