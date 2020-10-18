import numpy as np


def point_on_line_by_z(pnt1: np.ndarray, pnt2: np.ndarray, z: float)\
        -> np.ndarray:
    ratio = (z - pnt1[2]) / (pnt2[2] - pnt1[2])
    return pnt1 + (pnt2-pnt1)*ratio
