#include "scatter.h"

scatter::scatter(int z1, double m1, int z2, double m2)
{
    double m1_m2;

    m1_m2 = m1 / m2;
    rnorm = 0.4685 / (pow(z1, 0.23) + pow(z2, 0.23));
    enorm = 14.39979 * z1 * z2 / rnorm * (1 + m1_m2);
    dirfac = 2 / (1 + m1_m2);
    denfac = 4 * m1_m2 / ((1 + m1_m2) * (1 + m1_m2));
}

void scatter::ZBLscreen(double r, double* screen, double* dscreen)
{
    double exp1, exp2, exp3, exp4;

    exp1 = exp(-B1 * r);
    exp2 = exp(-B2 * r);
    exp3 = exp(-B3 * r);
    exp4 = exp(-B4 * r);

    *screen = A1*exp1 + A2*exp2 + A3*exp3 + A4*exp4;
    *dscreen = - (A1B1*exp1 + A2B2*exp2 + A3B3*exp3 + A4B4*exp4);
}

double scatter::estimate_apsis(double e, double p)
{
    double r0;
    double psq, r0sq;
    double screen, dscreen;
    double numerator, denominator;
    double residuum;

    psq = p*p;
    r0sq = 0.5 * (psq + sqrt(psq*psq + 4*K3/e));

    if(r0sq < R23sq){
        r0sq = psq + K2/e;
        if(r0sq < R12sq){
            r0 = (1 + sqrt(1 + 4*e*(e+K1)*psq)) / (2*(e+K1));
        } else{
            r0 = sqrt(r0sq);
        }
    } else{
        r0 = sqrt(r0sq);
    }

    // Do Newton-Raphson iterations to improve the estimate
    for(int i=0; i<NITER; i++){
        ZBLscreen(r0, &screen, &dscreen);
        numerator= r0 * (r0 - screen/e) - (p*p);
        denominator = 2*r0 - (screen+r0+dscreen)/e;
        r0 -= numerator/denominator;

        residuum = 1 - screen/(e*r0) - (p*p)/(r0*r0);
        if(abs(residuum) < 1e-4){
            break;
        }
    }

    return r0;
}

double scatter::magic(double e, double p)
{
    double r0, screen, dscreen;
    double rho, sqrte;
    double alpha, beta, gamma, delta, a, b, g;
    double cos_half_theta;

    r0 = estimate_apsis(e, p);
    ZBLscreen(r0, &screen, &dscreen);

    rho = 2*(e*r0-screen) / (screen/r0-dscreen);
    sqrte = sqrt(e);
    alpha = 1 + C1/sqrte;
    beta = (C2+sqrte) / (C3+sqrte);
    gamma = (C4+e) / (C5+e);
    a = 2 * alpha * e * pow(p, beta);
    g = gamma / (sqrt(1+(a*a))-a);
    delta = a * (r0-p) / (1+g);

    cos_half_theta = (p + rho + delta) / (r0 + rho);
    if(cos_half_theta > 1){
        std::cout << "Warning: cos_half_theta > 1: " << cos_half_theta << "\n";
        std:: cout << "e =" << e << "p =" << p << "r0 =" << r0 << "rho =" << rho << "delta =" << delta << "\n";
    }

    return cos_half_theta;
}

void scatter::doScatter(double* e, double dir[3], double p, double dirp[3], double dir_new[3], double dir_recoil[3], double* e_recoil)
{
    double cos_half_theta;
    double sin_psi, cos_psi;
    double norm = 0;

    // scattering angle theta in the center-of-mass system
    cos_half_theta = magic(*e/enorm, p/rnorm);

    // directions of the recoil and the projectile after the collision
    sin_psi = cos_half_theta;
    cos_psi = sqrt(1 - sin_psi*sin_psi);
    for(int i=0; i<3; i++)
    {
        dir_recoil[i] = dirfac * cos_psi * (cos_psi * dir[i] + sin_psi * dirp[i]);
        dir_new[i] = dir[i] - dir_recoil[i];
        
        norm += (dir_new[i] * dir_new[i]);
    }
    norm = sqrt(norm);

    if(norm == 0){
        std::copy(dir, dir+3, dir_new);
    } else {
        for(int i=0; i<3; i++){
            dir_new[i] = dir_new[i] / norm;
        }
    }

    norm = sqrt(dir_recoil[0]*dir_recoil[0] + dir_recoil[1]*dir_recoil[1] + dir_recoil[2]*dir_recoil[2]);
    if(norm == 0){
        std::copy(dir, dir+3, dir_recoil);
    } else{
        for(int i=0; i<3; i++){
            dir_recoil[i] = dir_recoil[i] / norm;
        }
    }

    *e_recoil = denfac * *e * (1 - cos_half_theta*cos_half_theta);
    *e -= *e_recoil;
}