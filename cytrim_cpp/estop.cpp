#include "estop.h"

estop::estop(double corr_lindhard, int z1, double m1, int z2, double density)
{
    this->density = density;

    this->fac_lindhard = corr_lindhard * 1.212 * pow(z1, 7/6.0) * z2 / (pow((pow(z1, 2/3.0) + pow(z2, 2/3.0)), 3/2.0) * sqrt(m1));
}

double estop::eloss(double e, double free_path)
{
    double dee;

    dee = fac_lindhard * density * sqrt(e) * free_path;

    if(dee > e){
        dee = e;
    }

    return dee;
}