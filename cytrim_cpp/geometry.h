#pragma once

class geometry{
    private:
        double zmin, zmax;

    public:
        geometry(double zmin, double zmax);
        bool is_inside_target(double pos[3]);
};