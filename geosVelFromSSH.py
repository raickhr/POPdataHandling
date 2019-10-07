from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from operators import getGrad
from readWrite import *


filename = 'flt_HIGH_PASS.pop.h.0009-01-05.nc'
gridFileName = 'flt_HIGH_PASS.pop.h.0009-01-05.nc'
fieldFileNameList = 'flt_HIGH_PASS.pop.h.0009-01-05.nc'

dataset = Dataset(filename)

OMEGA = np.array(dataset.variables['omega'])
g = np.array(dataset.variables['grav'])
rho = np.array(dataset.variables['rho_sw'])


dimensions, U_Grid_Var, T_Grid_Var = readGridInfo(gridFileName)

Zlen = dimensions[0]
Ylen = dimensions[1]
Xlen = dimensions[2]

var = T_Grid_Var.loc[T_Grid_Var['name'] == 'DXT']['val']
heading = var.keys()[0]
DX = var[heading]

var = T_Grid_Var.loc[T_Grid_Var['name'] == 'DYT']['val']
heading = var.keys()[0]
DY = var[heading]

var = T_Grid_Var.loc[T_Grid_Var['name'] == 'TLAT']['val']
heading = var.keys()[0]
phi = var[heading]

var = T_Grid_Var.loc[T_Grid_Var['name'] == 'KMT']['val']
heading = var.keys()[0]
KMT = var[heading]

def writeWithGeosVel(readFilename):

    fieldsDF, timeDict = readField(readFilename, ['SSH','RHO','SALT','TAUX','TAUY','time'])

    dateTime = timeDict['val']
    timeUnits = timeDict['units']
    timeCalendar = timeDict['calendar']

    print('read file at date \n{0:04d}-{1:02d}-{2:02d}'.format(
        dateTime.year, dateTime.month, dateTime.day))

    var = fieldsDF.loc[fieldsDF['name'] == 'SSH']['val']
    heading = var.keys()[0]
    SSH = var[heading]
    SSH = SSH[0,:,:]

    landMask = np.ma.getmask(np.ma.masked_where(KMT < 1, KMT))
    eqMask = np.ma.getmask(np.ma.masked_where(abs(phi) <= 3, phi))
    nanMask = np.ma.getmask(np.ma.masked_where(abs(SSH) > 1e10, SSH))

    mask = landMask + eqMask + nanMask

    SSH = np.ma.array(SSH,mask = mask,fill_value=float('nan')).filled()

    dSSH_dx, dSSH_dy = getGrad(SSH, DX, DY)

    u_gos = -g/(2*OMEGA*np.sin(phi)) * dSSH_dy
    v_gos = g/(2*OMEGA*np.sin(phi)) * dSSH_dx

    plt.pcolormesh(u_gos)
    plt.clim(-200,200)
    plt.show()

    xlen = np.shape(u_gos)[1]
    ylen = np.shape(u_gos)[0]

    ugos = np.empty((1,ylen,xlen))
    vgos = np.empty((1, ylen, xlen))

    ugos[0,:,:] = u_gos
    vgos[0,:,:] = v_gos

    u_gos_units = 'centimeter/sec'
    v_gos_units = 'centimeter/sec'

    u_gos_long_name = 'Zonal geostrophic velocity from SSH'
    v_gos_long_name = 'Meridional geostrophic velocity from SSH'

    ugos_Dict = {
        'name': 'UGOS',
        'long_name': u_gos_long_name,
        'units': u_gos_units,
        'val': ugos,
    }

    vgos_Dict = {
        'name': 'VGOS',
        'long_name': v_gos_long_name,
        'units': v_gos_units,
        'val': vgos,
    }

    appendDF = pd.DataFrame(data=[ugos_Dict,vgos_Dict])

    fieldsDF = fieldsDF.append(appendDF,ignore_index=True)

    writeNetcdf_withTime('GeostrophicVel_', Xlen, Ylen,
                         Zlen, dateTime, timeUnits, timeCalendar, fieldsDF)


writeWithGeosVel('flt_HIGH_PASS.pop.h.0009-01-05.nc')









# phi = np.radians(np.ma.array(phi, mask=eqMask,
#                              fill_value=np.float('nan')).filled())


# landMask = np.ma.getmask(np.ma.masked_where(KMT < 1, KMT))
# eqMask = np.ma.getmask(np.ma.masked_where(abs(phi) <= 3, phi))
# nanMask = np.ma.getmask(np.ma.masked_where(abs(SSH) > 1e10, SSH))

# SSH = np.array(dataset.variables['SSH'])[0, :, :]



# SSH = np.ma.array(SSH, mask=landMask + nanMask + eqMask,
#                   fill_value=np.float('nan')).filled()

# dSSH_dx, dSSH_dy = getGrad(SSH, DX, DY)

# u_gos = -g/(2*OMEGA*np.sin(phi)) * dSSH_dy
# v_gos = g/(2*OMEGA*np.sin(phi)) * dSSH_dx

