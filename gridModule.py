from netCDF4 import Dataset
import numpy as np
import pandas as pd

from readWrite import *

gridFileName = 'flt_HIGH_PASS.pop.h.0009-01-05.nc'

dimensions, U_Grid_Var, T_Grid_Var = readGridInfo(gridFileName)

Zlen = dimensions[0]
Ylen = dimensions[1]
Xlen = dimensions[2]

var = U_Grid_Var.loc[U_Grid_Var['name'] == 'DXU']['val']
heading = var.keys()[0]
DXU = var[heading]

var = U_Grid_Var.loc[U_Grid_Var['name'] == 'DYU']['val']
heading = var.keys()[0]
DYU = var[heading]

var = U_Grid_Var.loc[U_Grid_Var['name'] == 'ULAT']['val']
heading = var.keys()[0]
ULAT = var[heading]

var = T_Grid_Var.loc[T_Grid_Var['name'] == 'KMT']['val']
heading = var.keys()[0]
KMT = var[heading]
