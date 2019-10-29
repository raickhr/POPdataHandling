from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gc

from readWrite import *
from gridModule import *
from constants import *
from operators import getUAREA3dmid, getLandMaskT3d

averagedFile = glob(AVG_outpath+'POP_time_averaged.nc')[0]

landMaskT = np.ma.getmask(np.ma.masked_where(KMT < 1, KMT))
landMaskT3d = []#np.ones((Zlen,Ylen,Xlen), dtype =bool)

landMaskT3d = getLandMaskT3d()

avg_SALT_k0 = np.empty((1,Ylen,Xlen),dtype=float)

fieldsDF = readField(averagedFile, ['SALT','RHO'])

var = fieldsDF.loc[fieldsDF['name'] == 'SALT']['val']
heading = var.keys()[0]
avg_SALT = var[heading] / 1000  ## changing to msu units from psu
avg_SALT_k0[0, :, :] = avg_SALT[0, :, :]

var = fieldsDF.loc[fieldsDF['name'] == 'RHO']['val']
heading = var.keys()[0]
avg_RHO = var[heading]
avg_RHO = np.ma.array(avg_RHO, mask = landMaskT3d, fill_value=float('nan')).filled()

area_avg_RHO = np.empty((1, Zlen, Ylen, Xlen), dtype=float)
UAREA3D = getUAREA3dmid()
for k in range(Zlen):
    area_avg_RHO[0, k, :, :] = np.nansum(avg_RHO[k, :, :] * UAREA3D[k, :, :]) \
        / np.nansum(UAREA3D[k, :, :])

def writeWithProdFields(readFilename,writeFilePath):
    global avg_SALT_k0, area_avg_RHO
    fieldsDF, timeDict = readField(
        readFilename, ['UVEL',
                       'VVEL',
                       'WVEL',                       
                       'RHO',
                       'SHF',
                       'EVAP_F',
                       'PREC_F',
                       'TAUX', 
                       'TAUY', 
                       'time'])

    dateTime = timeDict['val']
    timeUnits = timeDict['units']
    timeCalendar = timeDict['calendar']

    print('Reading file at date \n{0:04d}-{1:02d}-{2:02d}'.format(
        dateTime.year, dateTime.month, dateTime.day))

    var = fieldsDF.loc[fieldsDF['name'] == 'UVEL']['val']
    heading = var.keys()[0]
    UVEL = var[heading]
    uvel_k0 = np.empty((1, Ylen, Xlen), dtype=float)
    uvel_k0[0,:,:] = UVEL[0,0,:,:]

    var = fieldsDF.loc[fieldsDF['name'] == 'VVEL']['val']
    heading = var.keys()[0]
    VVEL = var[heading]
    vvel_k0 = np.empty((1, Ylen, Xlen), dtype=float)
    vvel_k0[0, :, :] = VVEL[0, 0,:, :]

    var = fieldsDF.loc[fieldsDF['name'] == 'WVEL']['val']
    heading = var.keys()[0]
    WVEL = var[heading]

    var = fieldsDF.loc[fieldsDF['name'] == 'RHO']['val']
    heading = var.keys()[0]
    RHO = var[heading]

    rho_k0 = np.empty((1, Ylen, Xlen), dtype=float)
    rho_k0[0,:,:] = RHO[0,0,:,:]

    var = fieldsDF.loc[fieldsDF['name'] == 'SHF']['val']
    heading = var.keys()[0]
    SHF = var[heading] * 1e7/10000 ## changing from watts/m^2 to ergs/(sec cm^2)

    var = fieldsDF.loc[fieldsDF['name'] == 'PREC_F']['val']
    heading = var.keys()[0]
    PREC_F = var[heading] * 1000/10000 ## changing from kg/m^2 to gram/cm^2

    var = fieldsDF.loc[fieldsDF['name'] == 'EVAP_F']['val']
    heading = var.keys()[0]
    EVAP_F = var[heading] * 1000/10000  # changing from kg/m^2 to gram/cm^2

    var = fieldsDF.loc[fieldsDF['name'] == 'TAUX']['val']
    heading = var.keys()[0]
    TAUX = var[heading]

    var = fieldsDF.loc[fieldsDF['name'] == 'TAUY']['val']
    heading = var.keys()[0]
    TAUY = var[heading]

    RHO_UVEL = RHO * UVEL
    RHO_VVEL = RHO * VVEL
    RHO_WVEL = RHO * WVEL

    RHO_star = RHO - area_avg_RHO
    RHO_star_RHO_star = RHO_star * RHO_star

    UVEL_UVEL = UVEL*UVEL
    UVEL_VVEL = UVEL*VVEL
    UVEL_WVEL = UVEL*WVEL

    VVEL_VVEL = VVEL*VVEL
    VVEL_WVEL = VVEL*WVEL

    Js = SHF/(cp_sw * rho)
    Gs = avg_SALT_k0 * (EVAP_F - PREC_F)/rho_fw

    Js_RHO = Js * rho_k0
    Gs_RHO = Gs * rho_k0

    TAUX_UVEL = TAUX * uvel_k0
    TAUY_VVEL = TAUY * vvel_k0

    all_data_dict = []

    var_dict =  {
        'name': 'RHO_UVEL',
        'long_name': 'product of RHO and UVEL',
        'units': 'gram/(cm^2 sec)',
        'val': RHO_UVEL,
    }
    del RHO_UVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'RHO_VVEL',
        'long_name': 'product of RHO and VVEL',
        'units': 'gram/(cm^2 sec)',
        'val': RHO_VVEL,
    }

    del RHO_VVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'RHO_WVEL',
        'long_name': 'product of RHO and WVEL',
        'units': 'gram/(cm^2 sec)',
        'val': RHO_WVEL,
    }

    del RHO_WVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'RHO_star',
        'long_name': 'differece of RHO and areaAveraged RHO at the level',
        'units': 'gram/cm^3',
        'val': RHO_star,
    }

    del RHO_star
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'RHO_star_RHO_star',
        'long_name': 'differece of RHO and areaAveraged RHO at the level',
        'units': 'gram^2/cm^6',
        'val': RHO_star_RHO_star,
    }

    del RHO_star_RHO_star
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'UVEL_UVEL',
        'long_name': 'product of UVEL and UVEL',
        'units': 'cm^2/sec^2',
        'val': UVEL_UVEL,
    }


    del UVEL_UVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'UVEL_VVEL',
        'long_name': 'product of UVEL and VVEL',
        'units': 'cm^2/sec^2',
        'val': UVEL_VVEL,
    }

    del UVEL_VVEL
    gc.collect()
    all_data_dict.append(var_dict)


    var_dict = {
        'name': 'UVEL_WVEL',
        'long_name': 'product of UVEL and WVEL',
        'units': 'cm^2/sec^2',
        'val': UVEL_WVEL,
    }

    del UVEL_WVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'VVEL_VVEL',
        'long_name': 'product of VVEL and VVEL',
        'units': 'cm^2/sec^2',
        'val': VVEL_VVEL,
    }

    del VVEL_VVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'VVEL_WVEL',
        'long_name': 'product of VVEL and WVEL',
        'units': 'cm^2/sec^2',
        'val': VVEL_WVEL,
    }
    del VVEL_WVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'Js_RHO',
        'long_name': 'product of Js=SHF/(cp*rho_sw) and RHO',
        'units': '(Kelvin cm)/sec',
        'val': Js_RHO,
    }
    del Js_RHO
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'Gs_RHO',
        'long_name': 'product of Gs=avg_SALT_k0(EVAP-PREC)/rho_fw and RHO',
        'units': '(msu cm)/sec',
        'val': Gs_RHO,
    }
    del Gs_RHO
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'TAUX_UVEL',
        'long_name': 'product of TAUX and UVEL',
        'units': 'erg/(cm^2sec)',
        'val': TAUX_UVEL,
    }
    del TAUX_UVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'TAUY_VVEL',
        'long_name': 'product of TAUY and VVEL',
        'units': 'erg/(cm^2sec)',
        'val': TAUY_VVEL,
    }
    del TAUY_VVEL
    gc.collect()
    all_data_dict.append(var_dict)

    del UVEL,VVEL,WVEL,RHO,PREC_F,EVAP_F
    del TAUX,TAUY,Js,Gs,rho_k0,uvel_k0,vvel_k0
    gc.collect()

    writeDF = pd.DataFrame(data=all_data_dict)

    writeNetcdf_withTime(writeFilePath + 'ProductFields_', Xlen, Ylen,
                         Zlen, dateTime, timeUnits, timeCalendar, writeDF)
