from netCDF4 import Dataset
import numpy as np
import sys
from glob import glob

#POP_outpath = ("../POP2_output/")
POP_outpath = ("../TEST/")
PRD_outpath = ("../PRD_output/")
ENERGY_resrv = ("../ENG_output/")
AVG_outpath = ("../AVG_output/")

filename = glob(POP_outpath + '*')[0]

dataset = Dataset(filename)

OMEGA = np.array(dataset.variables['omega'])
g = np.array(dataset.variables['grav'])
rho = np.array(dataset.variables['rho_sw'])
rho_fw = np.array(dataset.variables['rho_fw'])
cp_sw = np.array(dataset.variables['cp_sw'])

alpha = -2.55e-4 # this is used from POP2 equation of State 
beta = 7.64e-1 # this is used from POP2 equation of State


fill_value = 9.96921e+36

