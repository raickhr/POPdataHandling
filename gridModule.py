from netCDF4 import Dataset
from constants import *
import numpy as np
import pandas as pd
from glob import glob

from readWrite import *

gridFileName = glob(POP_outpath+'*')[0]

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

var = U_Grid_Var.loc[U_Grid_Var['name'] == 'dz']['val']
heading = var.keys()[0]
dz = var[heading]

var = U_Grid_Var.loc[U_Grid_Var['name'] == 'ULAT']['val']
heading = var.keys()[0]
ULAT = var[heading]

var = U_Grid_Var.loc[U_Grid_Var['name'] == 'ULONG']['val']
heading = var.keys()[0]
ULONG = var[heading]

var = U_Grid_Var.loc[U_Grid_Var['name'] == 'UAREA']['val']
heading = var.keys()[0]
UAREA = var[heading]

var = U_Grid_Var.loc[U_Grid_Var['name'] == 'z_w_top']['val']
heading = var.keys()[0]
z_w_top = var[heading]

var = U_Grid_Var.loc[U_Grid_Var['name'] == 'z_w_bot']['val']
heading = var.keys()[0]
z_w_bot = var[heading]


var = T_Grid_Var.loc[T_Grid_Var['name'] == 'KMT']['val']
heading = var.keys()[0]
KMT = var[heading]

var = T_Grid_Var.loc[T_Grid_Var['name'] == 'TAREA']['val']
heading = var.keys()[0]
TAREA = var[heading]

var = T_Grid_Var.loc[T_Grid_Var['name'] == 'z_t']['val']
heading = var.keys()[0]
z_t = var[heading]
