import sys
import os
from glob import glob
from readWrite import *
from geosVelFromSSH import writeWithGeosVel
import pandas as pd


fileList = ['flt_HIGH_PASS.pop.h.0009-01-05.nc']

for fileName in fileList:
    outFileName = fileName.rstrip('.nc')+'_add.nc'
    cmd = 'ncap2 -O -s "UV = UVEL*VVEL" '+fileName+' '+outFileName





yearList = [1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006]
for year in yearList:
    filename1 = str(year)+'/extracted/'+str(year)+'_windASC.nc'
    filename2 = str(year)+'/extracted/'+str(year)+'_windDSC.nc'
    cmd1 = 'python fillBlank.py --var1=eastward_wind --var2=northward_wind --filename=' + \
        filename1
    cmd2 = 'python fillBlank.py --var1=eastward_wind --var2=northward_wind --filename=' + \
        filename2

    print(cmd1)
    os.system(cmd1)

    print(cmd2)
    os.system(cmd2)

### MAKE GEOSVEL FOR THE LIST OF FILES ###

filenameList = ['flt_HIGH_PASS.pop.h.0009-01-05.nc']

for inputFilename in filenameList:
    writeWithGeosVel(inputFilename)


### READ AVEARAGED FILE #####

averagedFile = ''

fieldsDF = readField( averagedFile, ['UGOS', 'VGOS', 'RHO', 'SALT', 'TAUX', 'TAUY'])

var = fieldsDF.loc[fieldsDF['name'] == 'UGOS']['val']
heading = var.keys()[0]
avgUGOS = var[heading]

var = fieldsDF.loc[fieldsDF['name'] == 'VGOS']['val']
heading = var.keys()[0]
avgVGOS = var[heading]

var = fieldsDF.loc[fieldsDF['name'] == 'SALT']['val']
heading = var.keys()[0]
avgSALT = var[heading]

var = fieldsDF.loc[fieldsDF['name'] == 'RHO']['val']
heading = var.keys()[0]
avgRHO = var[heading]


