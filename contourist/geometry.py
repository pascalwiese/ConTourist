import numpy as np
# import matplotlib.pyplot as plt


def point_on_line_by_z(pnt1: np.ndarray, pnt2: np.ndarray, z: float)\
        -> np.ndarray:
    ratio = (z - pnt1[2]) / (pnt2[2] - pnt1[2])
    pnt = pnt1 + (pnt2-pnt1)*ratio

    # plt.plot([pnt1[0], pnt2[0]], [pnt1[1], pnt2[1]], 'bo')
    # plt.scatter([pnt[0]], [pnt[1]], color='red')
    # plt.text(pnt1[0], pnt1[1], str(pnt1[2]))
    # plt.text(pnt2[0], pnt2[1], str(pnt2[2]))
    # plt.text(pnt[0], pnt[1], str(pnt[2]))
    # plt.show()

    return pnt
