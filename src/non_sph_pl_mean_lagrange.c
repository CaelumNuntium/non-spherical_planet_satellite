#include <math.h>
#include "NSP_mean_lagrange.h"

void mean_lagrange(double kappa_sq, double J2, double r_enc, double a, double e, double incl, double* M_dot, double* M0_dot, double* g_dot, double* Omega_dot)
{
    double mean_motion, s, c, eta, mu;
    mu = -J2;
    mean_motion = sqrt(kappa_sq) / (a * sqrt(a));
    eta = sqrt(1.0 - e * e);
    c = cos(incl);
    s = sin(incl);
    *M0_dot = mean_motion * mu * 0.75 * r_enc * r_enc / (a * a * eta * eta * eta) * (1.0 - 3.0 * c * c);
    *M_dot = mean_motion * (1.0 + mu * 0.75 * r_enc * r_enc / (a * a * eta * eta * eta) * (1.0 - 3.0 * c * c));
    *g_dot = mu * 0.75 * r_enc * r_enc * mean_motion / (a * a * eta * eta * eta * eta) * (1.0 - 5.0 * c * c);
    *Omega_dot = mu * mean_motion * 1.5 * r_enc * r_enc * c / (a * a * eta * eta * eta * eta);
}
