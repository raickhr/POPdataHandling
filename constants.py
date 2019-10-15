from netCDF4 import Dataset
import numpy as np

filename = 'flt_HIGH_PASS.pop.h.0009-01-05.nc'

dataset = Dataset(filename)

OMEGA = np.array(dataset.variables['omega'])
g = np.array(dataset.variables['grav'])
rho = np.array(dataset.variables['rho_sw'])
rho_fw = np.array(dataset.variables['rho_fw'])
cp_sw = np.array(dataset.variables['cp_sw'])
