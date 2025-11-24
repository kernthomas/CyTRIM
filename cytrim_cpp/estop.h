#pragma once

#include <iostream>
#include <cmath>

class estop{
    private:
        double fac_lindhard;
        double density;

    public:
        estop(double corr_lindhard, int z1, double m1, int z2, double density);
        double eloss(double e, double free_path);
};