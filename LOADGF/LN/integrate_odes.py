# *********************************************************************
# FUNCTION TO INTEGRATE ODEs for SPHEROIDAL DEFORMATION OF AN ELASTIC BODY
#
# *********************************************************************

import numpy as np
import math
import sys
from LOADGF.LN import integrate_fullEarth_n0
from LOADGF.LN import integrate_fullEarth
from LOADGF.LN import integrate_mantle
from LOADGF.LN import apply_boundary_conditions_n0
from LOADGF.LN import apply_boundary_conditions
from LOADGF.LN import evaluate_load_ln_n0
from LOADGF.LN import evaluate_load_ln
from LOADGF.LN import compute_solutions_n0
from LOADGF.LN import compute_solutions
from scipy import interpolate

def main(n,s_min,tck_lnd,tck_mnd,tck_rnd,tck_gnd,wnd,ond,piG,sic,soc,small,\
    num_soln,backend,abs_tol,rel_tol,nstps,order,gnd,adim,gsdim,L_sc,T_sc,inf_tol,s,nmaxfull):

    # Special Case: n=0
    # if (n == 0):
    #
    #     # # Perform the Integration Through Full Earth
    #     Y1, Y2, sint_mt = integrate_fullEarth_n0.main(n,s_min,tck_lnd,tck_mnd,tck_rnd,tck_gnd,wnd,ond,piG,sic,soc,small,\
    #         num_soln,backend,abs_tol,rel_tol,nstps,order,adim,gsdim,L_sc,T_sc,inf_tol,s)
    #     # # Apply Boundary Conditions at Surface
    #     m_load,m_pot,m_str,m_shr = apply_boundary_conditions_n0.main(n,Y1[-1],Y2[-1],gnd[-1],piG)
    #     # Compute Y Solutions
    #     Y_load = compute_solutions_n0.main(Y1,Y2,m_load)
    #     # Compute Load Love Numbers
    #     nkprime = evaluate_load_ln_n0.main(n)

    # else:
        # Perform the Integration Through the Mantle Only
    Y1, Y2, Y3, sint_mt = integrate_mantle.main(n,tck_lnd,tck_mnd,tck_rnd,tck_gnd,wnd,ond,piG,\
        num_soln,backend,abs_tol,rel_tol,nstps,order,inf_tol,s,soc)
    # Apply Boundary Conditions at Surface
    m_load,m_pot,m_str,m_shr = apply_boundary_conditions.main(n,Y1[-1],Y2[-1],Y3[-1],gnd[-1],piG)
    # Compute Y Solutions
    Y_load = compute_solutions.main(Y1,Y2,Y3,m_load)
    # Compute Load Love Numbers
    hprime,nlprime,nkprime = evaluate_load_ln.main(n,Y_load[-1,:],adim,gsdim,T_sc,L_sc)

    # Return Love Numbers
    nk = nkprime/n
    print('THIS IS nk: ', nk)

    return nkprime


