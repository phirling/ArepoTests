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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', type=int, default=2, choices=[1,2])
    args = parser.parse_args()

    # Main physical parameters
    BOXSIZE_PHYS_KPC = 1.0
    NDENS0_PHYS_CGS = 0.1
    N_STAR = 1
    STAR_MASS = 400
    GFM_N_CHEM_ELEMENTS = 10
    GFM_SOLAR_METALLICITY = 0.0127
    GLASS_IC_PATH = Path('/Users/patrickhirling/Galspec/misc/ICs/Glass_N16.hdf5')
    TEMPERATURE = 1e4 # K
    MAXNUMSNE=1000
    MU = 1.26 # (We assume an ionized primordial gas here, approx)
    GAMMA = 5.0 / 3.0
    MODE = int(args.mode) # 1: GFM/TNG, 2: MinGFM/SF_ECOGAL
    if MODE == 1:
        IC_OUTPUT_FILENAME = 'IC_GFMStar_N16.hdf5'
    else:
        IC_OUTPUT_FILENAME = 'IC_NewStar_N16.hdf5'

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
    BoxSize = FloatType(BOXSIZE_PHYS_KPC)
    mdens = FloatType(NDENS0_PHYS_CGS * m_p_cgs / UnitDensity_in_g_per_cm3)
    target_total_mass = BoxSize**3 * mdens # In internal mass units (Msun)

    # For cosmological run
    TimeBegin = 0.1

    # -------------------- Setup Gas
    # Load Glass Cube ICs
    with h5py.File(GLASS_IC_PATH) as f:
        Pos = np.array(f['PartType0/Coordinates'])
        Mass = np.array(f['PartType0/Masses'])
        Npart_gas = len(Pos)

    # Rescale the glass ICs to desired values
    Pos *= BOXSIZE_PHYS_KPC
    Mass *= target_total_mass

    # Set internal energy
    cell_u = ( k_b_cgs*TEMPERATURE / ((GAMMA-1) * MU * m_p_cgs) ) / UnitEnergy * UnitMass_in_g
    Utherm = np.full(Npart_gas, cell_u, dtype=FloatType)

    # Velocities
    Vel = np.zeros((Npart_gas, 3), dtype=FloatType)

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

    # -------------------- Setup Star
    # Create needed arrays
    Npart_stars = 1
    Pos_stars = np.zeros((Npart_stars, 3), dtype=FloatType)
    Vel_stars = np.zeros((Npart_stars, 3), dtype=FloatType)
    Mass_stars = np.zeros((Npart_stars), dtype=FloatType)
    stellar_formation_times = np.zeros(Npart_stars, dtype=np.float32)
    stellar_parent_densities = np.zeros(Npart_stars, dtype=np.float32)
    stellar_masses = np.zeros((Npart_stars, MAXNUMSNE), dtype=np.float32)
    stellar_explosion_times = np.zeros((Npart_stars, MAXNUMSNE), dtype=np.float32)
    stellar_NSNE = np.zeros(Npart_stars, dtype=IntType)

    # Enrichment
    stellar_hsml = np.zeros(Npart_stars, dtype=np.float32)
    stellar_let = np.zeros(Npart_stars, dtype=np.float32)
    stellar_Z = np.zeros(Npart_stars, dtype=np.float32)
    stellar_ZMass = np.zeros((Npart_stars, GFM_N_CHEM_ELEMENTS), dtype=np.float32)
    stellar_InitialMass = np.zeros((Npart_stars, GFM_N_CHEM_ELEMENTS), dtype=np.float32)

    # Explicitely set needed values
    Pos_stars[0,:] = BoxSize / 2.0

    Mass_stars[0] = FloatType(STAR_MASS)

    stellar_formation_times[0] = TimeBegin
    stellar_parent_densities[0] = mdens

    stellar_NSNE[0] = 1
    stellar_masses[0,0] = 8 / HubbleParam
    stellar_explosion_times[0,0] = 1.001 * TimeBegin

    stellar_hsml[0] = 1.5 * BoxSize / Npart_gas**(1./3) # ~1.5 x cell radius, VERY approximate first guess
    stellar_let[0] = TimeBegin
    stellar_InitialMass[0] = Mass_stars[0]

    # Here we use solar values for this example SSP, based on test_stellar_evolution in GFM
    stellar_Z[0]  = 1 * GFM_SOLAR_METALLICITY
    stellar_ZMass[0,0] = Mass_stars[0] * 0.76
    stellar_ZMass[0,1] = Mass_stars[0] * (1.0 - 0.76)
    stellar_ZMass[0,2] = stellar_Z[0] / GFM_SOLAR_METALLICITY * Mass_stars[0] * 2.06e-3
    stellar_ZMass[0,3] = stellar_Z[0] / GFM_SOLAR_METALLICITY * Mass_stars[0] * 8.36e-4
    stellar_ZMass[0,4] = stellar_Z[0] / GFM_SOLAR_METALLICITY * Mass_stars[0] * 5.49e-3
    stellar_ZMass[0,5] = stellar_Z[0] / GFM_SOLAR_METALLICITY * Mass_stars[0] * 1.41e-3
    stellar_ZMass[0,6] = stellar_Z[0] / GFM_SOLAR_METALLICITY * Mass_stars[0] * 5.91e-4
    stellar_ZMass[0,7] = stellar_Z[0] / GFM_SOLAR_METALLICITY * Mass_stars[0] * 6.83e-4
    stellar_ZMass[0,8] = stellar_Z[0] / GFM_SOLAR_METALLICITY * Mass_stars[0] * 1.10e-3

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
        NumPart = np.array([Npart_gas, 0, 0, 0, Npart_stars, 0], dtype=IntType)
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
        if MODE == 2:
            gp_gas.create_dataset('ElementAbundances', data = Abund)
            gp_gas.create_dataset('DusttoGasRatio', data = DustRatio)

        # PartType4 Dataset (Stars)
        if N_STAR > 0:
            gp_star.create_dataset('ParticleIDs', data = np.arange(Npart_gas + 1, Npart_gas + Npart_stars + 1))
            gp_star.create_dataset('Coordinates', data = Pos_stars / a * h)
            gp_star.create_dataset('Masses', data = Mass_stars / h)
            gp_star.create_dataset('Velocities', data = Vel_stars / a**(0.5))
            if MODE == 1:
                gp_star.create_dataset('GFM_StellarFormationTime', data = stellar_formation_times)
                gp_star.create_dataset('GFM_InitialMass', data = stellar_InitialMass / h)
                gp_star.create_dataset('GFM_Metallicity', data = stellar_Z)
                gp_star.create_dataset('GFM_Metals', data = stellar_ZMass / Mass_stars)
            else:
                gp_star.create_dataset('PID', data = np.arange(Npart_gas, Npart_gas + Npart_stars))
                gp_star.create_dataset('StellarFormationTime', data = stellar_formation_times)
                gp_star.create_dataset('InitialMass', data = stellar_InitialMass / h)
                gp_star.create_dataset('StellarMetallicity', data = stellar_Z)
                gp_star.create_dataset('StellarMetalMasses', data = stellar_ZMass / Mass_stars)
                gp_star.create_dataset('ParentDensity', data = stellar_parent_densities * a**3 / h**2)
                gp_star.create_dataset('StellarHsml', data = stellar_hsml * h / a)

                #gp_star.create_dataset('NumberOfSupernovae', data = stellar_NSNE)
                #gp_star.create_dataset('StellarMasses', data = stellar_masses)
                #gp_star.create_dataset('StellarExplosionTimes', data = stellar_explosion_times)

                #gp_star.create_dataset('LastEnrichTime', data = stellar_let)


    print(f"Estimated TargetGasMass: {target_total_mass/Npart_gas * h:.3g}. Mean effective cell mass: {(Mass*h).mean():.3g}")
    print("done.")