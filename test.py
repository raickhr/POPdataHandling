from netCDF4 import Dataset
import matplotlib.pyplot as plt
from glob import glob
from constants import *
from gridModule import *
from operators import *

averagedFile = glob(AVG_outpath+'POP_time_averaged.nc')[0]

landMaskT = np.ma.getmask(np.ma.masked_where(KMT < 1, KMT))
landMaskT3d = []  # np.ones((Zlen,Ylen,Xlen), dtype =bool)

landMaskT3d = getLandMaskT3d()

fieldsDF = readField(averagedFile, ['SALT', 'RHO'])

var = fieldsDF.loc[fieldsDF['name'] == 'SALT']['val']
heading = var.keys()[0]
avg_SALT = var[heading] / 1000  # changing to msu units from psu
avg_SALT = np.ma.array(avg_SALT, mask=landMaskT3d,
                       fill_value=float('nan')).filled()

avg_SALT_k0 = avg_SALT[0, :, :]

var = fieldsDF.loc[fieldsDF['name'] == 'RHO']['val']
heading = var.keys()[0]
avg_RHO = var[heading]
avg_RHO = np.ma.array(avg_RHO, mask=landMaskT3d,
                      fill_value=float('nan')).filled()


PRD_fileList = glob(PRD_outpath+'*')
PRD_fileList.sort()

POP_fileList = glob(POP_outpath+'*')
POP_fileList.sort()

POP_file = POP_fileList[0]
PRD_file = PRD_fileList[0]

POP_data = Dataset(POP_file)
PRD_data = Dataset(PRD_file)

SHF = np.array(POP_data.variables['SHF'])[0,:,:] * 1e7/10000
EVAP_F = np.array(POP_data.variables['EVAP_F'])[0,:,:] * 1000/10000
PREC_F = np.array(POP_data.variables['PREC_F'])[0, :, :] * 1000/10000
RHO_k0 = np.array(POP_data.variables['RHO'])[0,0,:,:]

Js = SHF/(cp_sw * rho)
Gs = avg_SALT_k0 * (EVAP_F - PREC_F)/rho_fw

Js_RHO = Js * RHO_k0
Gs_RHO = Gs * RHO_k0

Js_RHO2 = np.array(PRD_data.variables['Js_RHO'])[0,:,:]
Gs_RHO2 = np.array(PRD_data.variables['Gs_RHO'])[0,:,:]

Js_RHO = np.ma.array(Js_RHO, mask= landMaskT, fill_value= float('nan')).filled()
Gs_RHO = np.ma.array(Gs_RHO, mask=landMaskT, fill_value=float('nan')).filled()
Js_RHO2 = np.ma.array(Js_RHO2, mask=landMaskT, fill_value=float('nan')).filled()
Gs_RHO2 = np.ma.array(Gs_RHO2, mask=landMaskT, fill_value=float('nan')).filled()

print(np.shape(Js_RHO), np.shape(Js_RHO2), np.shape(Gs_RHO), np.shape(Gs_RHO2))

plt.subplot(2,2,1)
plt.pcolormesh(Js_RHO)
plt.colorbar()

plt.subplot(2, 2, 2)
plt.pcolormesh(Js_RHO2)
plt.colorbar()

plt.subplot(2, 2, 3)
plt.pcolormesh(Gs_RHO)
plt.clim(np.nanmin(Gs_RHO), np.nanmax(Gs_RHO))
plt.colorbar()

plt.subplot(2, 2, 4)
plt.pcolormesh(Gs_RHO2)
plt.clim(np.nanmin(Gs_RHO), np.nanmax(Gs_RHO))
plt.colorbar()

plt.show()



