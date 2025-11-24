#pragma once
#include <cmath>
#include <iostream>

#define A1  0.18175
#define A2  0.50986
#define A3  0.28022
#define A4  0.02817

#define B1  3.1998
#define B2  0.94229
#define B3  0.4029
#define B4  0.20162

#define A1B1    A1 * B1
#define A2B2    A2 * B2
#define A3B3    A3 * B3
#define A4B4    A4 * B4

// Constants for apsis estimation for the ZBL potential
#define K2  0.38            // factor of the 1/R part
#define K3  7.2             // factor of the 1/R^3 part
#define K1  1/(4 * K2)
#define R12sq   (2 * K2) * (2 * K2)
#define R23sq   K3 / K2
#define NITER   1           // number of Newton-Raphson iterations

#define C1  0.99229
#define C2  0.011615
#define C3  0.007122
#define C4  14.813
#define C5  9.3066

class scatter{
    private:
        double enorm, rnorm, dirfac, denfac;
    public:
        scatter(int z1, double m1, int z2, double m2);
        
        void ZBLscreen(double r, double* screen, double* dscreen);
        double estimate_apsis(double e, double p);
        double magic(double e, double p);
        void doScatter(double* e, double dir[3], double p, double dirp[3], double dir_new[3], double dir_recoil[3], double* e_recoil);
};