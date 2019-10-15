import numpy as np

def getGrad(array, DX, DY):
    ## this function calculates gradients 
    # on U grid based on the field on the T grid

    i_plus = np.roll(array, -1, axis=1)
    j_plus = np.roll(array, -1, axis=0)

    i_plus_j_plus = np.roll(array, -1, axis=1)

    gradX = 0.5 * (i_plus + i_plus_j_plus - array - j_plus)/DX
    gradY = 0.5 * (j_plus + i_plus_j_plus - array - i_plus)/DY

    return gradX, gradY
