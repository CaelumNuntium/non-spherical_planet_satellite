#include <math.h>
#include "non_spherical_planet.h"

extern double r_enc;
extern double J2;
extern double kappa_sq;

void non_spherical_planet_field(double t, double* x, double* res)
{
    double mu_sq, r_sq, r_cube, w_xy, w_z, ww;
    double mu;
    r_sq = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    r_cube = r_sq * sqrt(r_sq);
    mu = -J2;
    w_xy = -kappa_sq * (7.5 * r_enc * r_enc * mu * x[2] * x[2] + r_sq * r_sq - 1.5 * r_enc * r_enc * mu * r_sq) / (r_cube * r_sq * r_sq);
    w_z = -kappa_sq * (7.5 * r_enc * r_enc * mu * x[2] * x[2] + r_sq * r_sq - 4.5 * r_enc * r_enc * mu * r_sq) / (r_cube * r_sq * r_sq);
    res[0] = x[3];
    res[1] = x[4];
    res[2] = x[5];
    res[3] = w_xy * x[0];
    res[4] = w_xy * x[1];
    res[5] = w_z * x[2];
}
