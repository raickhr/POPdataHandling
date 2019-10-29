from netCDF4 import Dataset
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gc

from readWrite import *
from gridModule import *
from constants import *

landMaskT = np.ma.getmask(np.ma.masked_where(KMT < 1, KMT))
landMaskT3d = []  # np.ones((Zlen,Ylen,Xlen), dtype =bool)

for k in range(Zlen):
    lmask = np.ma.getmask(np.ma.masked_where(KMT < k+1, KMT))
    landMaskT3d.append(lmask)

landMaskT3d = np.array(landMaskT3d)


def averagePOP_output(readPath, writeFilePath):
    fileList = glob(readPath+'*')
    total_Files = len(fileList)

    avg_UVEL = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)
    avg_VVEL = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)
    avg_WVEL = np.zeros((1 ,Zlen, Ylen, Xlen), dtype=float)
    avg_RHO = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)
    avg_PD = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)
    avg_SHF = np.zeros((1, Ylen, Xlen), dtype=float)
    avg_EVAP_F = np.zeros((1, Ylen, Xlen), dtype=float)
    avg_PREC_F = np.zeros((1, Ylen, Xlen), dtype=float)
    avg_TAUX = np.zeros((1, Ylen, Xlen), dtype=float)
    avg_TAUY = np.zeros((1, Ylen, Xlen), dtype=float)
    avg_SALT = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)

    for readFilename in fileList:
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
                           'SALT',
                           'PD',
                           'time'])

        dateTime = timeDict['val']
        timeUnits = timeDict['units']
        timeCalendar = timeDict['calendar']

        print('Reading file at date \n{0:04d}-{1:02d}-{2:02d}'.format(
            dateTime.year, dateTime.month, dateTime.day))

        var = fieldsDF.loc[fieldsDF['name'] == 'UVEL']['val']
        heading = var.keys()[0]
        UVEL = var[heading]
        avg_UVEL +=  UVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'VVEL']['val']
        heading = var.keys()[0]
        VVEL = var[heading]
        avg_VVEL += VVEL
        
        var = fieldsDF.loc[fieldsDF['name'] == 'WVEL']['val']
        heading = var.keys()[0]
        WVEL = var[heading]
        avg_WVEL += WVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'RHO']['val']
        heading = var.keys()[0]
        RHO = var[heading]
        avg_RHO += RHO

        var = fieldsDF.loc[fieldsDF['name'] == 'PD']['val']
        heading = var.keys()[0]
        PD = var[heading]
        avg_PD += PD

        var = fieldsDF.loc[fieldsDF['name'] == 'SALT']['val']
        heading = var.keys()[0]
        SALT = var[heading]
        avg_SALT += SALT

        var = fieldsDF.loc[fieldsDF['name'] == 'SHF']['val']
        heading = var.keys()[0]
        SHF = var[heading]
        avg_SHF += SHF

        var = fieldsDF.loc[fieldsDF['name'] == 'PREC_F']['val']
        heading = var.keys()[0]
        PREC_F = var[heading]
        avg_PREC_F += PREC_F

        var = fieldsDF.loc[fieldsDF['name'] == 'EVAP_F']['val']
        heading = var.keys()[0]
        EVAP_F = var[heading]
        avg_EVAP_F += EVAP_F

        var = fieldsDF.loc[fieldsDF['name'] == 'TAUX']['val']
        heading = var.keys()[0]
        TAUX = var[heading]
        avg_TAUX += TAUX

        var = fieldsDF.loc[fieldsDF['name'] == 'TAUY']['val']
        heading = var.keys()[0]
        TAUY = var[heading]
        avg_TAUY += TAUY

    avg_UVEL /= total_Files
    avg_VVEL /= total_Files
    avg_WVEL /= total_Files
    avg_RHO /= total_Files
    avg_PD /= total_Files
    avg_SHF /= total_Files
    avg_EVAP_F /= total_Files
    avg_PREC_F /= total_Files
    avg_TAUX /= total_Files
    avg_TAUY /= total_Files
    avg_SALT /= total_Files

    all_data_dict = []

    var_dict = {
        'name': 'UVEL',
        'long_name': 'velocity along i direction',
        'units': 'cm/sec)',
        'val': avg_UVEL[0,:,:,:],
    }

    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'VVEL',
        'long_name': 'velocity along j direction',
        'units': 'cm/sec)',
        'val': avg_VVEL[0, :, :, :],
    }

    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'WVEL',
        'long_name': 'velocity along k direction(upwards)',
        'units': 'cm/sec)',
        'val': avg_WVEL[0, :, :, :],
    }

    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'RHO',
        'long_name': 'density',
        'units': 'gram/cm^3',
        'val': avg_RHO[0, :, :, :],
    }

    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'PD',
        'long_name': 'Potential Density Ref to Surface',
        'units': 'gram/cm^3',
        'val': avg_PD[0, :, :, :],
    }

    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'SALT',
        'long_name': 'salinity',
        'units': 'gm/kg',
        'val': avg_SALT[0, :, :, :],
    }
    
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'SHF',
        'long_name': 'Total Surface Heat Flux, Including SW',
        'units': 'watt/m^2',
        'val': avg_SHF[0, :, :],
    }

    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'EVAP_F',
        'long_name': 'Evaporation Flux from Coupler',
        'units': 'kg/m^2/s',
        'val': avg_EVAP_F[0, :, :],
    }

    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'PREC_F',
        'long_name': 'Precipitation Flux from Cpl (rain+snow)',
        'units': 'kg/m^2/s',
        'val': avg_PREC_F[0, :, :],
    }

    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'TAUX',
        'long_name': 'Windstress in grid-x direction',
        'units': 'dyne/cm^2',
        'val': avg_TAUX[0, :, :],
    }

    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'TAUY',
        'long_name': 'Windstress in grid-y direction',
        'units': 'dyne/cm^2',
        'val': avg_TAUY[0, :, :],
    }

    all_data_dict.append(var_dict)
    

    writeDF = pd.DataFrame(data=all_data_dict)

    writeNetcdf(writeFilePath + 'POP_time_averaged.nc', Xlen, Ylen, Zlen, writeDF)


def averagePRD_output(readPath, writeFilePath):
    fileList = glob(readPath+'*')
    total_Files = len(fileList)

    avg_RHO_UVEL = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)
    avg_RHO_VVEL = np.zeros((1, Zlen, Ylen, Xlen), dtype=float) 
    avg_RHO_WVEL = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)

    avg_RHO_star = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)

    avg_UVEL_UVEL = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)
    avg_UVEL_VVEL = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)
    avg_UVEL_WVEL = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)

    avg_VVEL_VVEL = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)
    avg_VVEL_WVEL = np.zeros((1, Zlen, Ylen, Xlen), dtype=float)

    avg_Js_RHO = np.zeros((1, Ylen, Xlen), dtype=float) 
    avg_Gs_RHO = np.zeros((1, Ylen, Xlen), dtype=float) 

    avg_TAUX_UVEL = np.zeros((1, Ylen, Xlen), dtype=float)
    avg_TAUY_VVEL = np.zeros((1, Ylen, Xlen), dtype=float)


    for readFilename in fileList:
        fieldsDF, timeDict = readField(
            readFilename, ['RHO_UVEL',
                           'RHO_VVEL',
                           'RHO_WVEL',
                           'RHO_star',
                           'UVEL_UVEL',
                           'UVEL_VVEL',
                           'UVEL_WVEL',
                           'VVEL_VVEL',
                           'VVEL_WVEL',
                           'Js_RHO',
                           'Gs_RHO',
                           'TAUX_UVEL',
                           'TAUY_VVEL',
                           'time'])

        dateTime = timeDict['val']
        timeUnits = timeDict['units']
        timeCalendar = timeDict['calendar']

        print('Reading file at date \n{0:04d}-{1:02d}-{2:02d}'.format(
            dateTime.year, dateTime.month, dateTime.day))

        var = fieldsDF.loc[fieldsDF['name'] == 'RHO_UVEL']['val']
        heading = var.keys()[0]
        RHO_UVEL = var[heading]
        avg_RHO_UVEL += RHO_UVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'RHO_VVEL']['val']
        heading = var.keys()[0]
        RHO_VVEL = var[heading]
        avg_RHO_VVEL += RHO_VVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'RHO_WVEL']['val']
        heading = var.keys()[0]
        RHO_WVEL = var[heading]
        avg_RHO_WVEL += RHO_WVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'RHO_star']['val']
        heading = var.keys()[0]
        RHO_star = var[heading]
        avg_RHO_star += RHO_star

        var = fieldsDF.loc[fieldsDF['name'] == 'UVEL_UVEL']['val']
        heading = var.keys()[0]
        UVEL_UVEL = var[heading]
        avg_UVEL_UVEL += UVEL_UVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'UVEL_VVEL']['val']
        heading = var.keys()[0]
        UVEL_VVEL = var[heading]
        avg_UVEL_VVEL += UVEL_VVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'UVEL_WVEL']['val']
        heading = var.keys()[0]
        UVEL_WVEL = var[heading]
        avg_UVEL_WVEL += UVEL_WVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'VVEL_VVEL']['val']
        heading = var.keys()[0]
        VVEL_VVEL = var[heading]
        avg_VVEL_VVEL += VVEL_VVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'VVEL_WVEL']['val']
        heading = var.keys()[0]
        VVEL_WVEL = var[heading]
        avg_VVEL_WVEL += VVEL_WVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'Js_RHO']['val']
        heading = var.keys()[0]
        Js_RHO = var[heading]
        avg_Js_RHO += Js_RHO

        var = fieldsDF.loc[fieldsDF['name'] == 'Gs_RHO']['val']
        heading = var.keys()[0]
        Gs_RHO = var[heading]
        avg_Gs_RHO += Gs_RHO

        var = fieldsDF.loc[fieldsDF['name'] == 'TAUX_UVEL']['val']
        heading = var.keys()[0]
        TAUX_UVEL = var[heading]
        avg_TAUX_UVEL += TAUX_UVEL

        var = fieldsDF.loc[fieldsDF['name'] == 'TAUY_VVEL']['val']
        heading = var.keys()[0]
        TAUY_VVEL = var[heading]
        avg_TAUY_VVEL += TAUY_VVEL

    avg_RHO_UVEL /= total_Files
    avg_RHO_VVEL /= total_Files
    avg_RHO_WVEL /= total_Files
    avg_RHO_star /= total_Files
    avg_UVEL_UVEL /= total_Files
    avg_UVEL_VVEL /= total_Files
    avg_UVEL_WVEL /= total_Files
    avg_VVEL_VVEL /= total_Files
    avg_VVEL_WVEL /= total_Files
    avg_Js_RHO /= total_Files
    avg_Gs_RHO /= total_Files
    avg_TAUX_UVEL /= total_Files
    avg_TAUY_VVEL /= total_Files

    all_data_dict = []

    var_dict = {
        'name': 'RHO_UVEL',
        'long_name': 'product of RHO and UVEL',
        'units': 'gram/(cm^2 sec)',
        'val': avg_RHO_UVEL[0, :, :, :],
    }
    del avg_RHO_UVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'RHO_VVEL',
        'long_name': 'product of RHO and VVEL',
        'units': 'gram/(cm^2 sec)',
        'val': avg_RHO_VVEL[0, :, :, :],
    }

    del avg_RHO_VVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'RHO_WVEL',
        'long_name': 'product of RHO and WVEL',
        'units': 'gram/(cm^2 sec)',
        'val': avg_RHO_WVEL[0, :, :, :],
    }

    del avg_RHO_WVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'RHO_star',
        'long_name': 'differece of RHO and areaAveraged RHO at the level',
        'units': 'gram/cm^3',
        'val': avg_RHO_star[0, :, :, :],
    }

    del avg_RHO_star
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'UVEL_UVEL',
        'long_name': 'product of UVEL and UVEL',
        'units': 'cm^2/sec^2',
        'val': avg_UVEL_UVEL[0, :, :, :],
    }

    del avg_UVEL_UVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'UVEL_VVEL',
        'long_name': 'product of UVEL and VVEL',
        'units': 'cm^2/sec^2',
        'val': avg_UVEL_VVEL[0, :, :, :],
    }

    del avg_UVEL_VVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'UVEL_WVEL',
        'long_name': 'product of UVEL and WVEL',
        'units': 'cm^2/sec^2',
        'val': avg_UVEL_WVEL[0, :, :, :],
    }

    del avg_UVEL_WVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'VVEL_VVEL',
        'long_name': 'product of VVEL and VVEL',
        'units': 'cm^2/sec^2',
        'val': avg_VVEL_VVEL[0, :, :, :],
    }

    del avg_VVEL_VVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'VVEL_WVEL',
        'long_name': 'product of VVEL and WVEL',
        'units': 'cm^2/sec^2',
        'val': avg_VVEL_WVEL[0, :, :, :],
    }
    del avg_VVEL_WVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'Js_RHO',
        'long_name': 'product of Js=SHF/(cp*rho_sw) and RHO',
        'units': '(Kelvin cm)/sec',
        'val': avg_Js_RHO[0, :, :],
    }
    del avg_Js_RHO
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'Gs_RHO',
        'long_name': 'product of Gs=avg_SALT_k0(EVAP-PREC)/rho_fw and RHO',
        'units': '(msu cm)/sec',
        'val': avg_Gs_RHO[0, :, :],
    }
    del avg_Gs_RHO
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'TAUX_UVEL',
        'long_name': 'product of TAUX and UVEL',
        'units': 'erg/(cm^2sec)',
        'val': avg_TAUX_UVEL[0, :, :],
    }
    del avg_TAUX_UVEL
    gc.collect()
    all_data_dict.append(var_dict)

    var_dict = {
        'name': 'TAUY_VVEL',
        'long_name': 'product of TAUY and VVEL',
        'units': 'erg/(cm^2sec)',
        'val': avg_TAUY_VVEL[0,:,:],
    }
    del avg_TAUY_VVEL
    gc.collect()

    all_data_dict.append(var_dict)

    writeDF = pd.DataFrame(data=all_data_dict)

    writeNetcdf(writeFilePath + 'PRD_time_averaged.nc',
                Xlen, Ylen, Zlen, writeDF)



