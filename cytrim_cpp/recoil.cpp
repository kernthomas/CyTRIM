#include "recoil.h"

recoil::recoil(double density) : gen(std::random_device{}()), dist(0.0, 1.0)
{
    mean_free_path = pow(density, -1/3.0);
    pmax = mean_free_path / sqrt(M_PI);
}

double recoil::rand()
{
    return dist(gen);
}

void recoil::get_recoil_position(const dirVector<double, 3>& pos, const dirVector<double, 3>& dir, double* free_path, double* p, dirVector<double, 3>& dirp, dirVector<double, 3>& pos_recoil)
{
    dirVector<double, 3> pos_collision;
    double fi, cos_fi, sin_fi;
    double cos_alpha, sin_alpha;
    double cos_phi, sin_phi;
    unsigned int k, i, j;
    double norm;

    *free_path = mean_free_path;

    pos_collision = pos + (dir * *free_path);

    *p = pmax * sqrt(rand());
    //*p = pmax * sqrt(0.5);

    // Azimuthal angle fi
    fi = 2 * M_PI * rand();
    //fi = 2 * M_PI * 0.5;
    cos_fi = cos(fi);
    sin_fi = sin(fi);

    // Convert direction vector to polar angles
    k = dir.vectabs().argmin();
    i = (k + 1) % 3;
    j = (i + 1) % 3;
    cos_alpha = dir[k];
    sin_alpha = sqrt(dir[i] * dir[i] + dir[j] * dir[j]);
    cos_phi = dir[i] / sin_alpha;
    sin_phi = dir[j] / sin_alpha;

    // direction vector from collision point to recoil
    dirp[i] = cos_fi * cos_alpha * cos_phi - sin_fi * sin_phi;
    dirp[j] = cos_fi * cos_alpha * sin_phi + sin_fi * cos_phi;
    dirp[k] = -cos_fi * sin_alpha;
    norm = dirp.norm();
    dirp = dirp / norm;
    
    pos_recoil = pos_collision + (dirp * *p);
}


void recoil::testRand()
{
    std::cout << rand() << "\n";
}
