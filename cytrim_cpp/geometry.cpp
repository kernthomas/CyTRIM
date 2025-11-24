#include "geometry.h"

geometry::geometry(double zmin, double zmax)
{
    this->zmax = zmax;
    this->zmin = zmin;
}

bool geometry::is_inside_target(double pos[3])
{
    if(zmin <= pos[2] && pos[2] <= zmax){
        return true;
    } else {
        return false;
    }
}