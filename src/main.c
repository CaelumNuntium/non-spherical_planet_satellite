#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "collo.h"
#include "legendre.h"
#include "rk4.h"
#include "non_spherical_planet.h"
#include "NSP_mean_lagrange.h"
#include "vector.h"
#include "orbits.h"

#define N_PARAMS 14

int write_result(const char*, int, double, double*, double, double, double, double, double, double, double, double, double);

double r_enc, J2, kappa_sq;

int main()
{
    int i, integrator, nt, nm, m, s;
    double h, a0, e0, i0, M0, g0, Omega0, dM, dM0, dg, d_Omega;
    vector3 r, v;
    double x0[6];
    double* rv_res;
    double* params_res;
    double* c;
    const char* parameters[N_PARAMS] = { "a0", "e0", "i0", "M0", "g0", "Omega0", "integrator", "delta_t", "n_steps", "m", "J2", "kappa_sq", "enclosing_radius", "s" };
    int is_required[N_PARAMS] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0 };
    char* values[N_PARAMS];
    for(i = 0; i < N_PARAMS; i++)
    {
        values[i] = (char*)malloc(100 * sizeof(char));
    }
    read_config("config.ini", N_PARAMS, parameters, values, is_required);
    sscanf(values[0], "%lf", &a0);
    sscanf(values[1], "%lf", &e0);
    sscanf(values[2], "%lf", &i0);
    i0 = radians(i0);
    sscanf(values[3], "%lf", &M0);
    M0 = radians(M0);
    sscanf(values[4], "%lf", &g0);
    g0 = radians(g0);
    sscanf(values[5], "%lf", &Omega0);
    Omega0 = radians(Omega0);
    if(!strcmp(values[6], "collo"))
    {
        integrator = 1;
    }
    else if(!strcmp(values[6], "rk4"))
    {
        integrator = 2;
    }
    else
    {
        fprintf(stderr, "Error: Invalid value: integrator = %s", values[3]);
        return 1;
    }
    sscanf(values[7], "%lf", &h);
    sscanf(values[8], "%d", &nt);
    if(!strcmp(values[9], "#ND"))
    {
        m = 1;
    }
    else
    {
        sscanf(values[9], "%d", &m);
    }
    sscanf(values[10], "%lf", &J2);
    sscanf(values[11], "%lf", &kappa_sq);
    sscanf(values[12], "%lf", &r_enc);
    if(integrator == 1)
    {
        if(!strcmp(values[13], "#ND"))
        {
            fprintf(stderr, "Error: Parameter s is required for collocation method");
            return 2;
        }
        else
        {
            sscanf(values[13], "%d", &s);
        }
    }
    for(i = 0; i < N_PARAMS; i++)
    {
        free(values[i]);
    }
    nm = nt / m;
    coords_from_orbital_elements(sqrt(kappa_sq), a0, e0, i0, Omega0, g0, M0, &r, &v);
    vector_to_array(r, x0);
    vector_to_array(v, x0 + 3);
    rv_res = (double*)malloc((nm + 1) * 6 * sizeof(double));
    switch(integrator)
    {
        case 2:
            runge_kutta(6, nt, m, h, non_spherical_planet_field, x0, rv_res);
            break;
        case 1:
            c = (double*)malloc(s * sizeof(double));
            lobatto(s, c);
            collo(6, nt, m, h, non_spherical_planet_field, x0, s, c, rv_res);
            free(c);
    }
    params_res = (double*)malloc((nm + 1) * 6 * sizeof(double));
    calculate_orbital_parameters(nm, 1, sqrt(kappa_sq), rv_res, params_res);
    free(rv_res);
    mean_lagrange(kappa_sq, J2, r_enc, a0, e0, i0, &dM, &dM0, &dg, &d_Omega);
    write_result("result.dat", nm, h * m, params_res, a0, e0, i0, M0, dM, g0, dg, Omega0, d_Omega);
    free(params_res);
    return 0;
}

int write_result(const char* filename, int nt, double h, double* res_params, double a0, double e0, double i0, double M0, double dM, double g0, double dg, double Omega0, double d_Omega)
{
    FILE* out;
    int i, j;
    double t, Omega, g, M;
    if(!(out = fopen(filename, "w")))
    {
        return 1;
    }
    fprintf(out, "# t    a    e    i    Omega    g    M0\n");
    for(i = 0; i < nt; i++)
    {
        t = h * i;
        Omega = Omega0 + d_Omega * t;
        Omega = to_interval(degrees(Omega), 0.0, 360.0);
        g = g0 + dg * t;
        g = to_interval(degrees(g), 0.0, 360.0);
        M = M0 + dM * t;
        M = to_interval(degrees(M), 0.0, 360.0);
        fprintf(out, "%.10lf ", t);
        for(j = 0; j < 2; j++)
        {
            fprintf(out, "%.10lf ", res_params[i * 6 + j]);
        }
        for(j = 2; j < 5; j++)
        {
            fprintf(out, "%.10lf ", degrees(res_params[i * 6 + j]));
        }
        fprintf(out, "%.10lf ", to_interval(degrees(res_params[i * 6 + 5]), 0.0, 360.0));
        fprintf(out, "%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", a0, e0, degrees(i0), Omega, g, M);
    }
    fclose(out);
}
