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

void recoil::get_recoil_position(double pos[3], double dir[3], double* free_path, double* p, double dirp[3], double pos_recoil[3])
{
    double pos_collision[3];
    double fi, cos_fi, sin_fi;
    double cos_alpha, sin_alpha;
    double cos_phi, sin_phi;
    unsigned int k, i, j;
    double norm;

    *free_path = mean_free_path;

    for(int idx=0; idx<3; idx++)
    {
        pos_collision[idx] = pos[idx] + *free_path * dir[idx];
    }

    *p = pmax * sqrt(rand());

    // Azimuthal angle fi
    fi = 2 * M_PI * rand();
    cos_fi = cos(fi);
    sin_fi = sin(fi);

    // Convert direction vector to polar angles
    auto it = std::min_element(dir, dir + 3, [](double a, double b) { return std::abs(a) < std::abs(b); });
    k = argmin_abs(dir);
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
    norm = sqrt(dirp[0]*dirp[0] + dirp[1]*dirp[1] + dirp[2]*dirp[2]);
    for(int idx=0; idx<3; idx++)
    {
        dirp[idx] /= norm;
        pos_recoil[idx] = pos_collision[idx] + *p * dirp[idx];
    }
}


void recoil::testRand()
{
    std::cout << rand() << "\n";
}

int recoil::argmin_abs(double dir[3])
{
    int k = 0;
    double min_val = std::abs(dir[0]);

    for (int i = 1; i < 3; ++i) {
        double v = std::abs(dir[i]);
        if (v < min_val) {
            min_val = v;
            k = i;
        }
    }

    return k;
}