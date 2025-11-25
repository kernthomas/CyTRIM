#include <iostream>
#include <chrono>
#include <cmath>
#include <omp.h>
#include "estop.h"
#include "recoil.h"
#include "geometry.h"
#include "dirVector.h"
#include "trajectory.h"


int main(int argc, char*argv[])
{
    unsigned int nion = 1000;   // number of projectiles to simulate, default value

    if (argc >= 2) {
        try {
            // argv[1] is the first argument after the program name
            nion = static_cast<unsigned int>(std::stoul(argv[1]));
        } catch (const std::exception& e) {
            std::cerr << "Invalid nion value: " << argv[1] << "\n";
            return 1;
        }
    } else {
        std::cout << "No nion specified, using default: " << nion << "\n";
        std::cout << "Usage: " << argv[0] << " <nion>\n";
    }

    std::cout << "Simulating " << nion << " projectiles...\n";

    auto start = std::chrono::high_resolution_clock::now();

    double zmin = 0.0;          // minimum z coordinate of the target (A)
    double zmax = 4000.0;        // maximum z coordinate of the target (A)
    int z1 = 5;                 // atomic number of projectile
    double m1 = 11.009;         // mass of projectile (amu)
    int z2 = 14;                // atomic number of target
    double m2 = 28.086;         // mass of target atom (amu)
    double density = 0.04994;   // target density (atoms/A^3)
    double corr_lindhard = 1.5;

    estop estopDev(corr_lindhard, z1, m1, z2, density);
    recoil recoilDev(density);
    scatter scatterDev(z1, m1, z2, m2);
    geometry geometryDev(zmin, zmax);
    trajectory trajectoryDev(recoilDev, estopDev, scatterDev, geometryDev);

    // Initial conditions of the projectile
    double e_init = 50000.0;                        // energy (ev)
    dirVector<double, 3> pos_init{0.0, 0.0, 0.0};   // position (A)
    dirVector<double, 3> dir_init{0.0, 0.0, 1.0};   // direction (unit vector)

    unsigned int count_inside = 0;
    double mean_z = 0;
    double std_z = 0;
    #pragma omp parallel
    {
        // Each thread runs its chunk of the loop
        #pragma omp for reduction(+:count_inside, mean_z, std_z)
        for(int i=0; i<nion; i++){
            double e;
            bool is_inside;
            dirVector<double, 3> pos;
            dirVector<double, 3> dir;

            trajectoryDev.simTrajectory(pos_init, dir_init, e_init, pos, dir, &e, &is_inside);

            if(is_inside){
                count_inside++;
                mean_z += pos[2];
                std_z += pos[2]*pos[2];
            }
        }
    }

    mean_z /= count_inside;
    std_z = sqrt(std_z / count_inside - mean_z*mean_z);

    std::cout << "Number of ions stopped inside the target: " << count_inside << "/" << nion << "\n";
    std::cout << "Mean penetration depth of ions stopped inside the target: " << mean_z << " A\n";
    std::cout << "Standard deviation of penetration depth: " << std_z << " A\n";

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Runtime: " << elapsed.count() << " seconds\n";

    return 0;
}