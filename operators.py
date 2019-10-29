import numpy as np
from constants import *

def getGrad(array, DX, DY):
    ## this function calculates gradients 
    # on U grid based on the field on the T grid
    # or vice versa

    i_plus = np.roll(array, -1, axis=1)
    j_plus = np.roll(array, -1, axis=0)

    i_plus_j_plus = np.roll(array, -1, axis=1)

    gradX = 0.5 * (i_plus + i_plus_j_plus - array - j_plus)/DX
    gradY = 0.5 * (j_plus + i_plus_j_plus - array - i_plus)/DY

    return gradX, gradY


def getGrad3D(ARRAY_IN, DX, DY, dz):
    from gridModule import Zlen, Ylen, Xlen
    landMaskT3d = getLandMaskT3d()

    array = np.ma.array(ARRAY_IN, mask = landMaskT3d, fill_value=float('nan')).filled()

    i_plus = np.roll(array, -1, axis=2)
    j_plus = np.roll(array, -1, axis=1)
    k_minus = np.roll(array, -1, axis=0)
    k_minus[Zlen-1,:,:] = 0.0

    DZ = np.empty((Zlen,Ylen,Xlen),dtype=float)
    for k in range(Zlen):
        DZ[k,:,:] = dz[k]

    i_plus_j_plus = np.roll(array, -1, axis=1)

    gradX = 0.5 * (i_plus + i_plus_j_plus - array - j_plus)/DX
    gradY = 0.5 * (j_plus + i_plus_j_plus - array - i_plus)/DY
    gradZ = (array - k_minus)/DZ

    return gradX, gradY, gradZ


def calc_stratification(PD_in):
    from gridModule import dz, Zlen, Ylen, Xlen

    landMaskT3d = getLandMaskT3d()
    PD = np.ma.array(PD_in, mask=landMaskT3d,fill_value=float('nan')).filled()

    DZ = np.empty((Zlen, Ylen, Xlen), dtype=float)
    area_avg_PD = np.empty((Zlen, Ylen, Xlen), dtype=float)

    AREA3D = getUAREA3dmid()

    for k in range(Zlen):
        DZ[k, :, :] = dz[k]
        area_avg_PD[k,:,:] = np.nansum(
            PD[k, :, :] * AREA3D[k, :, :])/np.nansum(AREA3D[k, :, :])

    k_plus = np.roll(area_avg_PD, -1,axis=0)
    k_plus[Zlen-1,:,:] = 0.0

    gradZ = (k_plus - area_avg_PD)/DZ

    n0 = np.empty((Zlen, Ylen, Xlen), dtype=float)
    for k in range(Zlen):
        n0[k, :, :] = gradZ[k]

    return n0

def getLandMaskT3d():
    from gridModule import KMT,Zlen
    ### CREATE LAND MASKS ###
    landMaskT3d = []  # np.ones((Zlen,Ylen,Xlen), dtype =bool)

    for k in range(Zlen):
        lmask = np.ma.getmask(np.ma.masked_where(KMT < k+1, KMT))
        landMaskT3d.append(lmask)

    landMaskT3d = np.array(landMaskT3d)

    return landMaskT3d



def getUAREA3dmid():
    from gridModule import UAREA, z_t, Zlen

    landMaskT3d = getLandMaskT3d()

    ### CALCULATE 3D AREA ###
    earth_rad = 6.371e8
    UAREA_3D = []
    KK = UAREA/(earth_rad**2)

    for k in range(Zlen):
        z = earth_rad - z_t[k]
        area = KK * z**2
        UAREA_3D.append(area)

    UAREA_3D = np.array(UAREA_3D)

    UAREA_3D = np.ma.array(UAREA_3D, mask=landMaskT3d,
                       fill_value=float('nan')).filled()

    return UAREA_3D

def getU_VOLUME():
    from gridModule import dz,Zlen,Ylen,Xlen
    DZ = np.empty((Zlen, Ylen, Xlen), dtype=float)

    for k in range(Zlen):
        DZ[k, :, :] = dz[k]

    UAREA_3D = getUAREA3dmid()
    U_VOL = UAREA_3D * DZ

    return U_VOL


def get_eddyPE_2_meanPE(eddy_RHO_UVEL, eddy_RHO_VVEL, avg_RHO, n0):
    from gridModule import Zlen, Ylen, Xlen, DXU, DYU

    h_GradX_rho = np.empty((Zlen, Ylen, Xlen), dtype=float)
    h_GradY_rho = np.empty((Zlen, Ylen, Xlen), dtype=float)
    g3D = np.empty((Zlen, Ylen, Xlen), dtype=float)

    for k in range(Zlen):
        h_GradX_rho[k,:,:], h_GradY_rho[k,:,:] =  \
            getGrad(avg_RHO[k, :, :], DXU, DYU)

        g3D[k, :, :] = g

    landMaskT3d = getLandMaskT3d()
        
    CONV = - g3D/n0 * (eddy_RHO_UVEL * h_GradX_rho + 
           eddy_RHO_VVEL * h_GradY_rho)

    CONV = np.ma.array(CONV,mask=landMaskT3d,fill_value=float('nan')).filled()

    VOLUME = getU_VOLUME()

    CONV = CONV * VOLUME

    totalCONV = np.nansum(CONV)

    return totalCONV


def get_EKE_2_MKE(eddy_UVEL_UVEL,
                  eddy_UVEL_VVEL,
                  eddy_UVEL_WVEL,
                  eddy_VVEL_VVEL,
                  eddy_VVEL_WVEL,
                  avg_UVEL,
                  avg_VVEL):
    from gridModule import Zlen, Ylen, Xlen, DXU, DYU, dz

    GradX_UVEL = np.empty((Zlen, Ylen, Xlen), dtype=float)
    GradY_UVEL = np.empty((Zlen, Ylen, Xlen), dtype=float)
    GradZ_UVEL = np.empty((Zlen, Ylen, Xlen), dtype=float)

    GradX_VVEL = np.empty((Zlen, Ylen, Xlen), dtype=float)
    GradY_VVEL = np.empty((Zlen, Ylen, Xlen), dtype=float)
    GradZ_VVEL = np.empty((Zlen, Ylen, Xlen), dtype=float)

    GradX_UVEL, GradY_UVEL, GradZ_UVEL =  \
        getGrad3D(avg_UVEL, DXU, DYU, dz)

    GradX_VVEL, GradY_VVEL, GradZ_VVEL =  \
        getGrad3D(avg_VVEL, DXU, DYU, dz)

    first_term = rho * eddy_UVEL_UVEL * GradX_UVEL
    second_term = rho * eddy_UVEL_VVEL * GradY_UVEL
    third_term = rho * eddy_UVEL_WVEL * GradZ_UVEL

    fourth_term = rho * eddy_UVEL_VVEL * GradX_VVEL
    fifth_term = rho * eddy_VVEL_VVEL * GradY_VVEL
    sixth_term = rho * eddy_VVEL_WVEL * GradZ_VVEL

    CONV = first_term + second_term + third_term + fourth_term + fifth_term + sixth_term

    VOLUME = getU_VOLUME()

    CONV = CONV * VOLUME

    landMaskT3d = getLandMaskT3d()

    CONV = np.ma.array(CONV, mask=landMaskT3d,
                       fill_value=float('nan')).filled()

    totalCONV = np.nansum(CONV)

    return totalCONV


def get_MAPE_2_MKE(avg_RHO, avg_W):
    CONV = - g * avg_RHO * avg_W

    VOLUME = getU_VOLUME()

    CONV = CONV * VOLUME

    landMaskT3d = getLandMaskT3d()

    CONV = np.ma.array(CONV, mask=landMaskT3d,
                       fill_value=float('nan')).filled()

    totalCONV = np.nansum(CONV)

    return totalCONV

def get_EAPE_2_EKE(eddy_RHO_WVEL):

    CONV = - g * eddy_RHO_WVEL
    VOLUME = getU_VOLUME()

    CONV = CONV * VOLUME

    landMaskT3d = getLandMaskT3d()

    CONV = np.ma.array(CONV, mask=landMaskT3d,
                       fill_value=float('nan')).filled()

    totalCONV = np.nansum(CONV)

    return totalCONV

def get_GENR_MAPE(avg_Js,avg_Gs,avg_RHO_star_k0,n0_k0):
    from gridModule import UAREA, KMT
    term1 = -g * alpha/n0_k0 * avg_Js * avg_RHO_star_k0 * UAREA
    term2 = -g * beta/n0_k0 * avg_Gs * avg_RHO_star_k0 * UAREA
    termSum = term1 + term2
    landMask = np.ma.getmask(np.ma.masked_where(KMT<1,KMT))
    termSum = np.ma.array(termSum, mask = landMask, fill_value=float('nan')).filled()
    GENR = np.nansum(termSum)
    return GENR


def get_GENR_EAPE(eddy_Js_RHO,eddy_Gs_RHO, n0_k0):
    from gridModule import UAREA, KMT
    

    term1 = -g * alpha/n0_k0 * eddy_Js_RHO * UAREA
    term2 = -g * beta/n0_k0 * eddy_Gs_RHO * UAREA
    termSum = term1 + term2
    landMask = np.ma.getmask(np.ma.masked_where(KMT < 1, KMT))
    termSum = np.ma.array(termSum, mask=landMask,
                          fill_value=float('nan')).filled()
    GENR = np.nansum(termSum)
    return GENR


def get_GENR_MKE(avg_TAUX,avg_TAUY,avg_UVEL,avg_VVEL):
    from gridModule import UAREA, KMT
    term1 = avg_TAUX * avg_UVEL * UAREA
    term2 = avg_TAUY * avg_VVEL * UAREA
    termSum = term1 + term2
    landMask = np.ma.getmask(np.ma.masked_where(KMT < 1, KMT))
    termSum = np.ma.array(termSum, mask=landMask,
                          fill_value=float('nan')).filled()
    GENR = np.nansum(termSum)
    return GENR


def get_GENR_EKE(eddy_TAUX_UVEL, eddy_TAUY_VVEL):
    from gridModule import UAREA, KMT
    termSum = (eddy_TAUX_UVEL + eddy_TAUY_VVEL) * UAREA
    landMask = np.ma.getmask(np.ma.masked_where(KMT < 1, KMT))
    termSum = np.ma.array(termSum, mask=landMask,
                          fill_value=float('nan')).filled()
    GENR = np.nansum(termSum)
    return GENR


def get_meanPE(avg_RHO_star,n0):
    from gridModule import Zlen, Ylen, Xlen
    g3D = np.empty((Zlen, Ylen, Xlen),dtype=float)
    for k in range(Zlen):
        g3D[k, :, :] = g

    ENG = g3D /(2*n0) * avg_RHO_star**2

    Volume = getU_VOLUME()
    landMaskT3d = getLandMaskT3d()

    ENG = ENG* Volume

    ENG = np.ma.array(ENG,mask = landMaskT3d, fill_value=float('nan')).filled()

    return np.nansum(ENG)



def get_meanKE(AVG_UVEL, AVG_VVEL):
    ENG = 1/2 * rho *(AVG_outpath**2 + AVG_VVEL**2)
    Volume = getU_VOLUME()
    landMaskT3d = getLandMaskT3d()

    ENG = ENG * Volume

    ENG = np.ma.array(ENG, mask=landMaskT3d, fill_value=float('nan')).filled()

    return np.nansum(ENG)

def get_eddyKE(eddy_UVEL_UVEL, eddy_VVEL_VVEL):
    ENG = 1/2 * rho * (eddy_UVEL_UVEL + eddy_VVEL_VVEL)
    Volume = getU_VOLUME()
    landMaskT3d = getLandMaskT3d()

    ENG = ENG * Volume

    ENG = np.ma.array(ENG, mask=landMaskT3d, fill_value=float('nan')).filled()

    return np.nansum(ENG)

def get_eddyPE(eddy_RHO_star,n0):
    from gridModule import Zlen, Ylen, Xlen
    g3D = np.empty((Zlen, Ylen, Xlen), dtype=float)
    for k in range(Zlen):
        g3D[k, :, :] = g

    ENG = g3D / (2*n0) * eddy_RHO_star

    Volume = getU_VOLUME()
    landMaskT3d = getLandMaskT3d()

    ENG = ENG * Volume

    ENG = np.ma.array(ENG, mask=landMaskT3d, fill_value=float('nan')).filled()

    return np.nansum(ENG)










    























