#pragma once
#include "dirVector.h"
#include "recoil.h"
#include "estop.h"
#include "scatter.h"
#include "geometry.h"

#define EMIN    5.0

class trajectory{
    private:
        recoil recoilDev;
        estop estopDev;
        scatter scatterDev;
        geometry geometryDev;

    public:
        trajectory(recoil& recoilDev, estop& estopDev, scatter& scatterDev, geometry& geometryDev):recoilDev(recoilDev), estopDev(estopDev), scatterDev(scatterDev), geometryDev(geometryDev){}

        void simTrajectory(const dirVector<double, 3>& pos_init, const dirVector<double, 3>& dir_init, const double e_init, dirVector<double, 3>& pos, dirVector<double, 3>& dir, double* e, bool* is_inside);

};