import numpy as np

def getGrad(array, DX, DY):
    i_minus = np.roll(array, 1, axis=1)
    i_plus = np.roll(array, -1, axis=1)

    j_minus = np.roll(array, 1, axis=0)
    j_plus = np.roll(array, -1, axis=0)

    two_dx = DX + 0.5 * np.roll(DX, 1, axis=1) + 0.5 * np.roll(DX, -1, axis=1)
    two_dy = DY + 0.5 * np.roll(DY, 1, axis=0) + 0.5 * np.roll(DY, -1, axis=0)

    gradX = (i_plus - i_minus)/(two_dx)
    gradY = (j_plus - j_minus)/(two_dy)

    return gradX, gradY
