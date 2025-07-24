import numpy as np
import astropy.units as au
import astropy.constants as ac
from astropy.cosmology import Planck15, z_at_value
import warnings
from pathlib import Path

"""
This is a useful script to determine parameters for cosmological tests, to have snapshots equally
spaced in proper time.
"""

def compute_time_params(z_init, SimTimeMyr, NumSnaps):
    """Compute the time/output parameters

    """
    SimTime     = float(SimTimeMyr) * au.Myr
    N_snaps     = int(NumSnaps)

    a_init = 1.0 / (z_init + 1.0)
    age_init = Planck15.age(z_init).to('Myr')
    age_end = age_init + SimTime
    z_end = float(z_at_value(Planck15.age, age_end))

    TimeBegin = a_init
    TimeMax = 1. / (z_end + 1)
    TimeBetSnapshot = (TimeMax / TimeBegin)**(1. / N_snaps)

    return TimeBegin, TimeMax, TimeBetSnapshot

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-z_init", type=float, default=9.0, help="Initial Redshift")
    parser.add_argument("-SimTimeMyr", type=float, default=100, help="Simulation time in proper Myr")
    parser.add_argument("-Nsnap", type=int, default=100, help="Number of snapshots")
    args = parser.parse_args()

    TimeBegin, TimeMax, TimeBetSnapshot = compute_time_params(args.z_init, args.SimTimeMyr, int(args.Nsnap))
    print(f"{'TimeBegin:':<30}{TimeBegin:.10g}")
    print(f"{'TimeMax:':<30}{TimeMax:.10g}")
    print(f"{'TimeBetSnapshot:':<30}{TimeBetSnapshot:.10g}")
