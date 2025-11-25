#pragma once
#include "dirVector.h"

class geometry{
    private:
        double zmin, zmax;

    public:
        geometry(double zmin, double zmax);
        bool is_inside_target(dirVector<double, 3>& pos);
};