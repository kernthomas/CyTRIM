#pragma once

#define _USE_MATH_DEFINES
#include <random>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "dirVector.h"

class recoil{
    private:
        double pmax;
        double mean_free_path;

        // Random number generator
        std::mt19937 gen;
        std::uniform_real_distribution<double> dist;
        double rand();

    public:
        recoil(double density);
        void get_recoil_position(const dirVector<double, 3>& pos, const dirVector<double, 3>& dir, double* free_path, double* p, dirVector<double, 3>& dirp, dirVector<double, 3>& pos_recoil);
        void testRand();
};