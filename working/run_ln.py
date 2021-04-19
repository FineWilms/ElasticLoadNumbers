#!/usr/bin/env python

# *********************************************************************
# MAIN PROGRAM TO COMPUTE LOVE NUMBERS (POTENTIAL/TIDE, LOAD, SHEAR, STRESS)
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

# IMPORT MPI MODULE
# from mpi4py import MPI

# MODIFY PYTHON PATH TO INCLUDE 'LoadDef' DIRECTORY
import sys
import os
sys.path.append(os.getcwd() + "/../")

# IMPORT PYTHON MODULES
import numpy as np
from LOADGF.LN import compute_love_numbers 

# --------------- SPECIFY USER INPUTS --------------------- #
 
# Full path to planet model text file
#     Planet model should be spherically symmetric, elastic, 
#         non-rotating, and isotropic (SNREI)
#     Format: radius(km), vp(km/s), vs(km/s), density(g/cc)
#     If the file delimiter is not whitespace, then specify in
#         call to function. 
planet_model = ("../input/Planet_Models/PREM.txt")
 
# Extension for the output filename (Default is '.txt')
file_ext      = (".txt")

# ------------------ END USER INPUTS ----------------------- #

# -------------------- BEGIN CODE -------------------------- #

# Ensure that the Output Directories Exist
# if (rank == 0):
if not (os.path.isdir("../output/Love_Numbers/")):
    os.makedirs("../output/Love_Numbers/")

# Compute Love Numbers
nk = compute_love_numbers.main(planet_model,file_out='wilms.txt')

# --------------------- END CODE --------------------------- #

