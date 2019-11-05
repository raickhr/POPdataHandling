from netCDF4 import Dataset
from netCDF4 import date2num, num2date
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gc

from readWrite import *
from gridModule import *
from constants import *
from operators import getUAREA3dmid, getLandMaskT3d

averagedFile = glob(AVG_outpath+'POP_time_averaged.nc')[0]

landMaskT = np.empty((1,Ylen,Xlen),dtype=bool)
landMaskT[0,:,:] = np.ma.getmask(np.ma.masked_where(KMT < 1, KMT))
landMaskT3d = []#np.ones((Zlen,Ylen,Xlen), dtype =bool)

landMaskT3d = getLandMaskT3d()

fieldsDF = readField(averagedFile, ['SALT','RHO'])

var = fieldsDF.loc[fieldsDF['name'] == 'SALT']['val']
heading = var.keys()[0]
avg_SALT = var[heading] / 1000  ## changing to msu units from psu
avg_SALT = np.ma.array(avg_SALT, mask=landMaskT3d,
                      fill_value=float('nan')).filled()

avg_SALT_k0 = np.empty((1, Ylen, Xlen), dtype=float)
avg_SALT_k0[0, :, :] = avg_SALT[0, :, :]

var = fieldsDF.loc[fieldsDF['name'] == 'RHO']['val']
heading = var.keys()[0]
avg_RHO = var[heading]
avg_RHO = np.ma.array(avg_RHO, mask = landMaskT3d, fill_value=float('nan')).filled()

area_avg_avg_RHO = np.empty((1, Zlen, Ylen, Xlen), dtype=float)
UAREA3D = getUAREA3dmid()

for k in range(Zlen):
    area_avg_avg_RHO[0, k, :, :] = np.nansum(avg_RHO[k, :, :] * UAREA3D[k, :, :]) \
        / np.nansum(UAREA3D[k, :, :])

def writeWithProdFields(readFilename,writeFilePath):
    global avg_SALT_k0, area_avg_RHO, landMaskT, landMaskT3d

    dataset = Dataset(readFilename)

    UVEL = np.array(dataset.variables['UVEL'])
    VVEL = np.array(dataset.variables['VVEL'])
    WVEL = np.array(dataset.variables['WVEL'])
    RHO = np.array(dataset.variables['RHO'])
    SHF = np.array(dataset.variables['SHF'])
    EVAP_F = np.array(dataset.variables['EVAP_F'])
    PREC_F = np.array(dataset.variables['PREC_F'])
    TAUX = np.array(dataset.variables['TAUX'])
    TAUY = np.array(dataset.variables['TAUY'])

    timeLongName = dataset.variables['time'].long_name
    timeUnits = dataset.variables['time'].units
    timeCalendar = dataset.variables['time'].calendar
    timeNum = np.array(dataset.variables['time'])[0]
    dateTime = num2date(timeNum, units=timeUnits, calendar=timeCalendar)

    print('Reading file at date \n{0:04d}-{1:02d}-{2:02d}'.format(
        dateTime.year, dateTime.month, dateTime.day))

    UVEL = np.ma.array(UVEL, mask=landMaskT3d,
                           fill_value=float('nan')).filled()
    VVEL = np.ma.array(VVEL, mask=landMaskT3d,
                           fill_value=float('nan')).filled()
    WVEL = np.ma.array(WVEL, mask=landMaskT3d,
                           fill_value=float('nan')).filled()
    RHO = np.ma.array(RHO, mask=landMaskT3d,
                          fill_value=float('nan')).filled()

    SHF = np.ma.array(SHF, mask=landMaskT,
                          fill_value=float('nan')).filled()

    PREC_F = np.ma.array(PREC_F, mask=landMaskT,
                         fill_value=float('nan')).filled()
    
    EVAP_F = np.ma.array(EVAP_F, mask=landMaskT,
                         fill_value=float('nan')).filled()
    TAUX = np.ma.array(TAUX, mask=landMaskT,
                           fill_value=float('nan')).filled()
    TAUY = np.ma.array(TAUY, mask=landMaskT,
                           fill_value=float('nan')).filled()

    uvel_k0 = np.empty((1, Ylen, Xlen), dtype=float)
    uvel_k0[0, :, :] = UVEL[0, 0, :, :]

    vvel_k0 = np.empty((1, Ylen, Xlen), dtype=float)
    vvel_k0[0, :, :] = VVEL[0, 0, :, :]

    rho_k0 = np.empty((1, Ylen, Xlen), dtype=float)
    rho_k0[0, :, :] = RHO[0, 0, :, :]

    RHO_UVEL = RHO * UVEL
    RHO_VVEL = RHO * VVEL
    RHO_WVEL = RHO * WVEL

    RHO_star = RHO - area_avg_avg_RHO
    RHO_star_RHO_star = RHO_star * RHO_star

    UVEL_UVEL = UVEL*UVEL
    UVEL_VVEL = UVEL*VVEL
    UVEL_WVEL = UVEL*WVEL

    VVEL_VVEL = VVEL*VVEL
    VVEL_WVEL = VVEL*WVEL

    SHF = SHF * 1e7/10000  # changing from watts/m^2 to ergs/(sec cm^2)
    EVAP_F = EVAP_F * 1000/10000  # changing from kg/m^2 to gram/cm^2
    PREC_F = PREC_F * 1000/10000  # changing from kg/m^2 to gram/cm^2

    Js = SHF/(cp_sw * rho)
    Gs = avg_SALT_k0 * (EVAP_F - PREC_F)/rho_fw

    Js_RHO = Js * rho_k0
    Gs_RHO = Gs * rho_k0

    TAUX_UVEL = TAUX * uvel_k0
    TAUY_VVEL = TAUY * vvel_k0

    RHO_UVEL = np.ma.array(
        RHO_UVEL, mask=landMaskT3d, fill_value=fill_value).filled()
    RHO_VVEL = np.ma.array(
        RHO_VVEL, mask=landMaskT3d, fill_value=fill_value).filled()
    RHO_WVEL = np.ma.array(
        RHO_WVEL, mask=landMaskT3d, fill_value=fill_value).filled()

    RHO_star = np.ma.array(
        RHO_star, mask=landMaskT3d, fill_value=fill_value).filled()
    RHO_star_RHO_star = np.ma.array(
        RHO_star_RHO_star, mask=landMaskT3d, fill_value=fill_value).filled()

    UVEL_UVEL = np.ma.array(
        UVEL_UVEL, mask=landMaskT3d, fill_value=fill_value).filled()
    UVEL_VVEL = np.ma.array(
        UVEL_VVEL, mask=landMaskT3d, fill_value=fill_value).filled()
    UVEL_WVEL = np.ma.array(
        UVEL_WVEL, mask=landMaskT3d, fill_value=fill_value).filled()
    VVEL_VVEL = np.ma.array(
        VVEL_VVEL, mask=landMaskT3d, fill_value=fill_value).filled()
    VVEL_WVEL = np.ma.array(
        VVEL_WVEL, mask=landMaskT3d, fill_value=fill_value).filled()

    Js_RHO = np.ma.array(Js_RHO, mask=landMaskT,
                             fill_value=fill_value).filled()
    Gs_RHO = np.ma.array(Gs_RHO, mask=landMaskT,
                             fill_value=fill_value).filled()

    TAUX_UVEL = np.ma.array(
        TAUX_UVEL, mask=landMaskT, fill_value=fill_value).filled()
    TAUY_VVEL = np.ma.array(
        TAUY_VVEL, mask=landMaskT, fill_value=fill_value).filled()


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


#writeWithProdFields('../TEST/rho_add_hist_test.pop.h.0009-01-05.nc',PRD_outpath)
