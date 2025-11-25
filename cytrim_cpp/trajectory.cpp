#include "trajectory.h"

void trajectory::simTrajectory(const dirVector<double, 3>& pos_init, const dirVector<double, 3>& dir_init, const double e_init, dirVector<double, 3>& pos, dirVector<double, 3>& dir, double* e, bool* is_inside)
{
    double freepath = 0;
    double p = 0;
    double dee;
    double e_recoil = 0;
    dirVector<double, 3> dirp;
    dirVector<double, 3> pos_recoil;
    dirVector<double, 3> dir_new;
    dirVector<double, 3> dir_recoil;

    pos = pos_init;
    dir = dir_init;
    *e = e_init;
    *is_inside = true;

    while(*e > EMIN){
        recoilDev.get_recoil_position(pos, dir, &freepath, &p, dirp, pos_recoil);
        dee = estopDev.eloss(*e, freepath);
        *e -= dee;

        pos = pos + (dir * freepath);
        if(geometryDev.is_inside_target(pos) == false){
            *is_inside = false;
            break;
        }

        scatterDev.doScatter(e, dir, p, dirp, dir_new, dir_recoil, &e_recoil);
        dir = dir_new;
    }
}