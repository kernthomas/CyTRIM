#include <iostream>
#include <cmath>
#include "estop.h"


int main()
{
    int z1 = 5;                     // atomic number of projectile
    double m1 = 11.009;         // mass of projectile (amu)
    int z2 = 14;                // atomic number of target
    double m2 = 28.086;         // mass of target atom (amu)
    double density = 0.04994;   // target density (atoms/A^3)
    double corr_lindhard = 1.5;

    estop estopInst(corr_lindhard, z1, m1, z2, density);

    return 0;
}