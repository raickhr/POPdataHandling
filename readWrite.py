from netCDF4 import Dataset
import numpy as np
import datetime
from netCDF4 import date2num, num2date
import pandas as pd
import sys

def readGridInfo(filename):
    dataset = Dataset(filename)

    dimensions = [dataset.dimensions['z_t'].size,
                  dataset.dimensions['nlat'].size,
                  dataset.dimensions['nlon'].size]

    U_GridVarNames = ['ULAT',
                      'ULONG',
                      'DXU',
                      'DYU',
                      'HUS',
                      'HUW',
                      'ANGLE',
                      'KMU']

    T_GridVarNames = ['TLAT',
                      'TLONG',
                      'DXT',
                      'DYT',
                      'HTE',
                      'HTN',
                      'ANGLET',
                      'KMT']

    UGridVar = []
    TGridVar = []

    for i in range(len(U_GridVarNames)):
        varName = U_GridVarNames[i]
        varLongName = dataset.variables[varName].long_name
        varUnits = ''
        if varName != 'KMU':
            varUnits = dataset.variables[varName].units
        val = np.array(dataset.variables[varName])

        varDict = {
            'name': varName,
            'long_name': varLongName,
            'units': varUnits,
            'val': val,
        }

        UGridVar.append(varDict)

    for i in range(len(T_GridVarNames)):
        varName = T_GridVarNames[i]
        varLongName = dataset.variables[varName].long_name
        varUnits = ''
        if varName != 'KMT':
            varUnits = dataset.variables[varName].units
        val = np.array(dataset.variables[varName])

        varDict = {
            'name': varName,
            'long_name': varLongName,
            'units': varUnits,
            'val': val,
        }

        TGridVar.append(varDict)

    

    U_Grid_Df = pd.DataFrame(data=UGridVar)
    T_Grid_Df = pd.DataFrame(data=TGridVar)

    return dimensions, U_Grid_Df, T_Grid_Df


def readField(filename, varNameList):
    dataset = Dataset(filename)

    readVars = []

    isTimeRequested = False

    for i in range(len(varNameList)):
        timeDict = {}
        varName = varNameList[i]
        varLongName = dataset.variables[varName].long_name
        varUnits = dataset.variables[varName].units
        val = np.array(dataset.variables[varName])

        if varName == 'time':
            isTimeRequested = True
            timeLongName = dataset.variables[varName].long_name
            timeUnits = dataset.variables[varName].units
            calendar = dataset.variables[varName].calendar
            timeNum = np.array(dataset.variables[varName])[0]
            time = num2date(timeNum,units=timeUnits,calendar=calendar)
            
            timeDict = {
                'name': 'time',
                'long_name': timeLongName,
                'units': timeUnits,
                'calendar': calendar,
                'val': time,
            }

        else:

            varDict = {
                'name': varName,
                'long_name': varLongName,
                'units': varUnits,
                'val': val,
            }

            readVars.append(varDict)

    fieldsDF = pd.DataFrame(data=readVars)

    if isTimeRequested:
        return fieldsDF, timeDict
    
    else:
        return fieldsDF


def writeNetcdf_withTime(writeFilename, Xlen, Ylen, Zlen, time, timeUnits, calendar, writeDF):

    year = time.year
    month = time.month
    day = time.day

    suffix = '{0:04d}-{1:02d}-{2:02d}.nc'.format(year,month,day)
    writeFilename = writeFilename + suffix

    writeDataset = Dataset(writeFilename, 'w', format='NETCDF4_CLASSIC')

    writeDataset.createDimension('nlon', Xlen)
    writeDataset.createDimension('nlat', Ylen)
    writeDataset.createDimension('z_t', Zlen)
    writeDataset.createDimension('time', None)

    timeVar = writeDataset.createVariable('time', np.float64, ('time'))
    timeVar.units = timeUnits

    timeVar[0] = date2num(time,timeUnits,calendar=calendar)

    numOfFields = len(writeDF)

    #writeVar = []

    for i in range(numOfFields):
        var = writeDF.iloc[[i]]['name']
        heading = var.keys()[0]
        name = var[heading]

        var = writeDF.iloc[[i]]['long_name']
        heading = var.keys()[0]
        long_name = var[heading]

        var = writeDF.iloc[[i]]['units']
        heading = var.keys()[0]
        units = var[heading]

        var = writeDF.iloc[[i]]['val']
        heading = var.keys()[0]
        array = var[heading]
        shape = np.shape(array)

        if len(shape) == 3:
            var = writeDataset.createVariable(name, float, ('time','nlat', 'nlon'))
            var.long_name = long_name
            var.units = units
            var[:, :, :] = array

        elif len(shape) == 4:
            var = writeDataset.createVariable(name, float, ('time','z_t','nlat', 'nlon'))
            var.long_name = long_name
            var.units = units
            var[:, :, :, :] = array

        else:
            print('{0} array size mismatch!!  size is {1}'.format(name, shape))
            sys.exit()


    writeDataset.close()


def writeNetcdf2D(writeFilename, Xlen, Ylen, time, timeUnits, writeDictList):
    writeDataset = Dataset(writeFilename, 'w', format='NETCDF4_CLASSIC')

    writeDataset.createDimension('nlon', Xlen)
    writeDataset.createDimension('nlat', Ylen)
    writeDataset.createDimension('time', None)

    numOfFields = len(writeDictList)

    #writeVar = []

    for i in range(numOfFields):
        name = writeDictList[i]['name']
        long_name = writeDictList[i]['long_name']
        units = writeDictList[i]['units']
        array = writeDictList[i]['val']

        var = writeDataset.createVariable(name, float, ('nlat', 'nlon'))
        var.long_name = long_name
        var.units = units
        var[:, :] = array

    writeDataset.close()
