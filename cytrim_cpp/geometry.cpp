#include "geometry.h"

geometry::geometry(double zmin, double zmax)
{
    this->zmax = zmax;
    this->zmin = zmin;
}

bool geometry::is_inside_target(dirVector<double, 3>& pos)
{
    if(zmin <= pos[2] && pos[2] <= zmax){
        return true;
    } else {
        return false;
    }
}