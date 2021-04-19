# *********************************************************************
# FUNCTION TO COMPUTE LOVE NUMBERS
#
# Copyright (c) 2014-2019: HILARY R. MARTENS, LUIS RIVERA, MARK SIMONS         
#
# This file is part of LoadDef.
#
#    LoadDef is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    LoadDef is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with LoadDef.  If not, see <https://www.gnu.org/licenses/>.
#
# *********************************************************************

# Import Python Modules
from __future__ import print_function
from mpi4py import MPI
import numpy as np
import math
import os
import sys
import datetime
# import matplotlib.pyplot as plt
from scipy import interpolate

# Import Modules from LoadDef
from LOADGF.LN import prepare_planet_model
# from LOADGF.LN import compute_asymptotic_LLN
from LOADGF.LN import integrate_odes
import pandas as pd
import matplotlib.pyplot as plt

"""
Compute load-deformation coefficients (also known as load Love numbers) based on an input 
spherically symmetric, non-rotating, elastic, and isotropic (SNREI) planetary model.
 
Input planetary model should be in the format [radius (km), Vp (km/s), Vs (km/s), density (g/cc)]

Parameters
----------
startn : Spherical harmonic degrees (for Love numbers) will be computed starting from this value
    Default is 0

stopn : Spherical harmonic degrees (for Love numbers) will be computed ending at this value
    Default is 10000

period_hours : Tidal forcing period (in hours)
    Default is 12.42 (M2 period)

r_min : Minimum radius for variable planetary structural properties (meters)
    Default is 1000

interp_emod : Optionally interpolate the planetary model to a different resolution
    Default is False (see LOADGF/LN/prepare_planet_model.py)

kx : Order of the spline fit for the planetary model (1=linear; 3=cubic)
    Default is 1 (recommended)

delim : Delimiter for the planetary model file
    Default is None (Whitespace)

inf_tol : Defines the integration starting radius, r, for which 
    :: (r/a)^n drops below influence tolerance, 'inf_tol'
    Default is 1E-5

rel_tol : Integration tolerance level (relative tolerance)
    Default is 1E-13

abs_tol : Integration tolerance level (absolute tolerance)
    Default is 1E-13

backend : Specify ODE Solver
    :: Recommended to only use 'dop853' or 'dopri5' solvers, since they integrate only to a specified stopping point (no overshoot, like lsoda and vode)
    Default is 'dop853'

nstps : Specify Maximum Number of (Internally Defined) Steps Allowed During Each Call to the ODE Solver
    :: For More Information, See Scipy.Integrate.Ode Manual Pages
    Default is 3000

num_soln : Set Number of Solutions for Each Integration (Integer)
    :: Note that integration step size is adaptive based on the
    ::  specified tolerance, but solutions are only computed at 
    ::  regular intervals determined by this user-specified value
    Default is 100

G : Universal Gravitational Constant
    Default is 6.672E-11 m^3/(kg*s^2)

nmaxfull : Maximum spherical harmonic degree for which integration will be performed through the full planet
           Beyond nmaxfull, integration will begin in the mantle
    Default is None (estimated from inf_tol within integrate_odes.py)

file_out : Extension for the output files.
    Default is ".txt"
"""


# Main Function
def main(myfile, startn=1, stopn=20, delim=None, period_hours=12.42, r_min=1000, inf_tol=1E-5,
         rel_tol=1E-12, abs_tol=1E-15, backend='dopri5', nstps=300, G=6.672E-11, file_out='wilms_nk.txt', kx=1, num_soln=1000,
         interp_emod=False, nmaxfull=1):
    # Print Status
    print(" ")
    print(":: Computing Love Numbers. Please Wait...")
    print(" ")
    # For SNREI Planet, Angular Frequency (omega) is Zero
    # Azimuthal order is only utilized for a rotating planet
    # The variables for each are included here as "place-holders" for future versions
    omega = 0
    order = 4

    # Prepare the planetary Model (read in, non-dimensionalize elastic parameters, etc.)
    r, mu, K, lmda, rho, g, tck_lnd, tck_mnd, tck_rnd, tck_gnd, s, lnd, mnd, rnd, gnd, s_min, small, \
    planet_radius, planet_mass, sic, soc, adim, gsdim, pi, piG, L_sc, R_sc, T_sc = \
        prepare_planet_model.main(myfile, G=G, r_min=r_min, kx=kx, file_delim=delim, emod_interp=interp_emod)
    # Define Forcing Period
    w = (1. / (period_hours * 3600.)) * (2. * pi)  # convert period to frequency (rad/sec)
    wnd = w * T_sc  # non-dimensionalize
    ond = omega * T_sc

    myn = np.linspace(startn, stopn, ((stopn - startn) + 1))
    n_sub = myn

    NK = []
    for ii in range(0,len(n_sub)):
        current_n = n_sub[ii]
        print('Working on Harmonic Degree: %7s | Number: %6d of %6d ' % (str(int(current_n)), (ii + 1), len(n_sub)))
        # Compute Integration Results for the Current Spherical Harmonic Degree, n

        nkprime = integrate_odes.main(current_n, s_min, tck_lnd, tck_mnd, tck_rnd, tck_gnd, wnd, ond, piG, sic, soc, small,
                            num_soln, backend, abs_tol, \
                            rel_tol, nstps, order, gnd, adim, gsdim, L_sc, T_sc, inf_tol, s, nmaxfull)


        NK.append(nkprime[0]/current_n)

    df = pd.DataFrame()
    df['nk'] = NK
    print(df)
    df.plot()
    plt.show()



    return
