# *********************************************************************
# FUNCTION TO COMPUTE LOAD LOVE NUMBERS FROM Y-SOLUTIONS (n>=1)
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

import numpy as np

def main(n,Y_load,a,gs,T_sc,L_sc):

    # Extract Solutions
    Y1sol_load = Y_load[0]
    Y2sol_load = Y_load[1]
    Y3sol_load = Y_load[2]
    Y4sol_load = Y_load[3]
    Y5sol_load = Y_load[4]
    Y6sol_load = Y_load[5]

    # Compute Load Love Numbers
    hprime = Y1sol_load
    nlprime = n*Y3sol_load
    Y5sol = Y5sol_load*((L_sc**2.)*(T_sc**(-2.)))
    nkprime = n*((Y5sol/(a*gs))-1.)

    # Adjust degree-one Love numbers to ensure that the potential field 
    # outside the Earth vanishes in the CE frame (e.g. Merriam 1985)
    if (n == 1):
        hprime = hprime - nkprime
        nlprime = nlprime - nkprime
        nkprime = nkprime - nkprime

    # Flatten Arrays
    hprime  = hprime.flatten()
    nlprime = nlprime.flatten()
    nkprime = nkprime.flatten()

    # Return Variables
    return hprime,nlprime,nkprime


