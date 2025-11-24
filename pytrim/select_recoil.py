"""Create the recoil position for the next collision.

Currently, only amorphous targets are supported. The free path length to
the next collision is assumed to be constant and equal to the atomic
density to the power -1/3.

Available functions:
    setup: setup module variables.
    get_recoil_position: get the recoil position.
"""
from math import sqrt, sin, cos, pi
import numpy as np
import cython


def setup(density):
    """Setup module variables depending on target density.

    Parameters:
        density (float): target density (atoms/A^3)

    Returns:
        None    
    """
    global PMAX, MEAN_FREE_PATH

    MEAN_FREE_PATH = density**(-1/3)
    PMAX = MEAN_FREE_PATH / sqrt(np.pi)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.locals(
    free_path=cython.double,
    p=cython.double,
    fi=cython.double,
    cos_fi=cython.double, sin_fi=cython.double,
    cos_alpha=cython.double, sin_alpha=cython.double,
    cos_phi=cython.double,  sin_phi=cython.double,
    norm=cython.double,
    i=cython.int, j=cython.int, k=cython.int,
    t=cython.int,
    di=cython.double, dj=cython.double,
    a0=cython.double, a1=cython.double, a2=cython.double,
)
@cython.ccall  # Python-callable when compiled; no-op under plain Python
def get_recoil_position(pos: np.ndarray, dir: np.ndarray):
    """
    Get the recoil position based on the projectile position and direction.

    Parameters
    ----------
    pos : ndarray, shape (3,)
    dir : ndarray, shape (3,)

    Returns
    -------
    free_path : float
    p         : float
    dirp      : ndarray, shape (3,)
    pos_recoil: ndarray, shape (3,)
    """
    if pos.shape[0] != 3 or dir.shape[0] != 3:
        raise ValueError("pos and dir must be length-3 arrays")

    free_path = MEAN_FREE_PATH

    # Random draws
    p  = PMAX * sqrt(float(np.random.random()))
    fi = 2.0 * pi * float(np.random.random())
    cos_fi = cos(fi)
    sin_fi = sin(fi)

    # Choose k = argmin(|dir|) so sin(alpha) stays large
    a0 = abs(dir[0])
    a1 = abs(dir[1])
    a2 = abs(dir[2])
    if a1 < a0:
        k = 1; a0 = a1
    else:
        k = 0
    if a2 < a0:
        k = 2

    i = (k + 1) % 3
    j = (i + 1) % 3

    cos_alpha = float(dir[k])
    di = float(dir[i]); dj = float(dir[j])
    sin_alpha = sqrt(di*di + dj*dj)

    # Guard against exactly axial directions
    if sin_alpha != 0.0:
        cos_phi = di / sin_alpha
        sin_phi = dj / sin_alpha
    else:
        cos_phi = 1.0
        sin_phi = 0.0

    # Recoil direction (unnormalized), then normalize
    dirp = np.empty(3, dtype=np.float64)
    dirp[i] =  cos_fi * cos_alpha * cos_phi - sin_fi * sin_phi
    dirp[j] =  cos_fi * cos_alpha * sin_phi + sin_fi * cos_phi
    dirp[k] = -cos_fi * sin_alpha

    norm = sqrt(dirp[0]*dirp[0] + dirp[1]*dirp[1] + dirp[2]*dirp[2])
    if norm != 0.0:
        dirp[0] /= norm; dirp[1] /= norm; dirp[2] /= norm

    # pos_recoil = (pos + free_path*dir) + p*dirp  (no broadcasting/temps)
    pos_recoil = np.empty(3, dtype=np.float64)
    for t in range(3):
        pos_recoil[t] = pos[t] + free_path * dir[t] + p * dirp[t]

    return free_path, p, dirp, pos_recoil
