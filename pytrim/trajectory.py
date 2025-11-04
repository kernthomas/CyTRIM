"""Simulate projectile trajectories.

Available functions:
    setup: setup module variables.
    trajectory: simulate one trajectory."""
from select_recoil import get_recoil_position
from scatter import scatter
from estop import eloss
from geometry import is_inside_target
import numpy as np
import cython

def setup():
    """Setup module variables.

    Parameters:
        None

    Returns:
        None    
    """
    global EMIN

    EMIN = 5.0  # eV


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.locals(
    e=cython.double,
    dee=cython.double,
    free_path=cython.double,
    is_inside=cython.bint,
    i=cython.int
)
@cython.ccall  # Python-callable when compiled; no-op under plain Python
def trajectory(pos_init: np.ndarray, dir_init: np.ndarray, e_init: cython.double):
    """Simulate one trajectory.
    
    Parameters:
        pos_init (ndarray): initial position of the projectile (size 3)
        dir_init (ndarray): initial direction of the projectile (size 3)
        e_init (float): initial energy of the projectile (eV)

    Returns:
        ndarray: final position of the projectile (size 3)
        ndarray: final direction of the projectile (size 3)
        float: final energy of the projectile (eV)
        bool: True if projectile is stopped inside the target, 
            False otherwise
    """
    pos = pos_init.copy()
    direction = dir_init.copy()
    e = e_init
    is_inside = True

    while e > EMIN:
        free_path, p, dirp, _ = get_recoil_position(pos, direction)

        dee = eloss(e, free_path)
        e -= dee

        # pos += free_path * direction   (no NumPy broadcasting/temporaries)
        for i in range(3):
            pos[i] += free_path * direction[i]

        if not is_inside_target(pos):
            is_inside = False
            break

        direction, e, _, _ = scatter(e, direction, p, dirp)

    return pos, direction, e, is_inside