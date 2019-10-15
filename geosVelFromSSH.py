from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from operators import getGrad
from readWrite import *
from gridModule import *
from constants import OMEGA, g

landMask = np.ma.getmask(np.ma.masked_where(KMT < 1, KMT))
eqMask = np.ma.getmask(np.ma.masked_where(abs(ULAT) <= 3, ULAT))

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

    nanMask = np.ma.getmask(np.ma.masked_where(abs(SSH) > 1e10, SSH))

    mask = landMask + eqMask + nanMask

    SSH = np.ma.array(SSH,mask = mask,fill_value=float('nan')).filled()
    phi = np.ma.array(ULAT, mask=mask, fill_value=float('nan')).filled()

    dSSH_dx, dSSH_dy = getGrad(SSH, DXU, DYU)

    u_gos = -g/(2*OMEGA*np.sin(np.radians(phi))) * dSSH_dy
    v_gos = g/(2*OMEGA*np.sin(np.radians(phi))) * dSSH_dx

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

