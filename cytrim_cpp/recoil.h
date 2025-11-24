#pragma once

#define _USE_MATH_DEFINES
#include <random>
#include <iostream>
#include <cmath>
#include <algorithm>

class recoil{
    private:
        double pmax;
        double mean_free_path;

        // Random number generator
        std::mt19937_64 gen;
        std::uniform_real_distribution<double> dist;
        double rand();

        int argmin_abs(double dir[3]);

    public:
        recoil(double density);
        void get_recoil_position(double* pos, double* dir, double* free_path, double* p, double* dirp, double* pos_recoil);
        void testRand();
};