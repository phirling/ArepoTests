import numpy as np
import h5py
import astropy.units as au
import argparse
from pathlib import Path
import warnings

def replace_multiple_lines(in_file_path, out_file_path, target_params, param_values):
    """Replace a list of parameters in txt file with given values
    """
    assert(isinstance(target_params,list))

    lines = []
    nparams = len(target_params)
    found_mask = np.zeros(nparams, dtype='bool')

    with open(in_file_path, 'r') as f:
        for i, line in enumerate(f):
            found_in_line = 0
            for k, param in enumerate(target_params):
                if param in line:
                    if found_mask[k]:
                        warnings.warn(f"Parameter '{target_params[k]}' was found twice in file")
                    else:
                        found_mask[k] = True
                    found_in_line += 1
                    if found_in_line > 1:
                        warnings.warn(f"Found more than one parameter on line {i:n}")

                    if type(param_values[k]) == str:
                        new_line = f"{param:<50}{param_values[k]}"
                    else:
                        new_line = f"{param:<50}{param_values[k]:.10g}"
                    lines.append(new_line + '\n')
            
            if found_in_line == 0:
                lines.append(line)

    if np.any(np.logical_not(found_mask)):
        warnings.warn("One or more parameters not found in file")

    with open(out_file_path, 'w') as f:
        f.writelines(lines)

# Main physical parameters
N_CELLS_PER_DIM = 64
BOXSIZE = 1.0
EPS = 1. / 100
fR = 1. / 4
TEMP = 1e4
gamma = 5.0 / 3.0
FloatType = np.float64
IntType = np.int32
mu = 1.26 # Primordial gas

# -------------------- Constants
m_p_cgs = 1.67262158e-24 # g
k_b_cgs = 1.380649e-16
G_cgs = 6.6743e-08 # cm3 / (g s2)
Msun_cgs = 1.9891e33

# -------------------- Units
UnitLength_in_cm            = 3.08568e21    # 1 kpc
UnitMass_in_g               = Msun_cgs     # 1 Solar Mass
UnitVelocity_in_cm_per_s    = 100000        # km / s
UnitTime_in_s               = UnitLength_in_cm / UnitVelocity_in_cm_per_s
UnitDensity_in_g_per_cm3 = UnitMass_in_g / UnitLength_in_cm**3
UnitEnergy = UnitMass_in_g * UnitLength_in_cm**2 / UnitTime_in_s**2
uts = UnitTime_in_s * au.s

# -------------------- Set Main Parameters
IC_filename = f"IC_UniformCloud_N{N_CELLS_PER_DIM:n}.hdf5"
param_filename = "param.txt"
Boxsize = FloatType(BOXSIZE)
N = IntType(N_CELLS_PER_DIM)

if __name__ == '__main__':

    Npart_gas = N**3

    Redshift = FloatType(0.0)
    Time = FloatType(0.0)

    dL = Boxsize / FloatType(N)
    dV = dL**3

    fact1 = (1 - EPS)*4./3 * np.pi * fR**3 + EPS
    invfact1 = 1.0 / fact1
    fact2 = 1 + EPS*(1 + 3.0/(4*np.pi*fR**3))
    invfact2 = 1.0 / fact2


    mdens_c = invfact1
    M_c = invfact2
    N_c = 4.0 / 3 * np.pi * fR**3 * Npart_gas

    mdens_a = EPS * mdens_c
    M_a = 1.0 - M_c

    # Estimate of target mass resolution
    cell_mass_target = M_c / N_c

    # Convert to absolute Path
    IC_filename = Path.cwd() / IC_filename
    param_filename = Path.cwd() / param_filename

    # ============================= Setup Gas ============================= #
    rng = np.random.default_rng(seed=1906) # For reproducibility
    Pos = rng.uniform(0.0, Boxsize, size=(Npart_gas, 3)).astype(FloatType)
    Rad = np.linalg.norm(Pos-Boxsize/2.0, axis=1)
    idx_incloud = np.where(Rad <= fR)[0]

    # Velocities
    Vel = np.zeros((Npart_gas, 3), dtype=FloatType)

    # "Masses" (In reality density)
    Mass = np.full(Npart_gas, mdens_a, dtype=FloatType)
    Mass[idx_incloud] = mdens_c

    # Internal Energy
    cell_u = ( k_b_cgs*TEMP / ((gamma-1) * mu * m_p_cgs) ) / UnitEnergy * UnitMass_in_g
    Utherm = np.full(Npart_gas, cell_u, dtype=FloatType)

    print(f"UnitTime                = {UnitTime_in_s:.4e} s  = {uts.to('Gyr'):.4e}")

    # ============================= Create Output File ============================= #
    print("Writing file: ", IC_filename, "...")
    with h5py.File(IC_filename, 'w') as f:
        header = f.create_group('Header')
        gp_gas = f.create_group('PartType0')
        gp_star = f.create_group('PartType4')

        # Header
        NumPart = np.array([Npart_gas, 0, 0, 0, 0, 0], dtype=IntType)
        header.attrs.create('NumPart_ThisFile', NumPart)
        header.attrs.create('NumPart_Total', NumPart)
        header.attrs.create('NumPart_Total_HighWord', np.zeros(6, dtype=IntType))
        header.attrs.create('MassTable', np.zeros(6, dtype=IntType))
        header.attrs.create('Time', Time)
        header.attrs.create('Redshift', Redshift)
        header.attrs.create('BoxSize', Boxsize)
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
        gp_gas.create_dataset('Coordinates', data = Pos)
        gp_gas.create_dataset('Masses', data = Mass)
        gp_gas.create_dataset('Velocities', data = Vel)
        gp_gas.create_dataset('InternalEnergy', data = Utherm)

    # ============================= Parameter File ============================= #
    params = ['ReferenceGasPartMass', 'InitCondFile']
    param_vals = [cell_mass_target, str(IC_filename.stem)]
    for i, pr in enumerate(params):
        pass
        #print(f"Editing file: {param_filename}: {pr} --> {param_vals[i]:.5g} ...")
    replace_multiple_lines(param_filename, param_filename, params, param_vals)