#include "vector.h"

void coords_from_orbital_elements(double kappa, double a, double e, double i, double Omega, double omega, double M, vector3* r, vector3* v);
void orbital_elements(double kappa, vector3 r, vector3 v, double* a, double* e, double* i, double* Omega, double* omega, double* M);
void calculate_orbital_parameters(int nt, int n_bodies, double kappa, double* coords, double* parameters);
double degrees(double);
double radians(double);
double to_interval(double value, double left_bound, double right_bound);
