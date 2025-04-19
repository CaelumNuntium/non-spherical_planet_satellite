#include <math.h>
#include "vector.h"
#include "orbits.h"

void coords_from_orbital_elements(double kappa, double a, double e, double i, double Omega, double omega, double M, vector3* r, vector3* v)
{
	const double eps = 2.5e-16;
	const int NMAX = 1000;
	int k;
	double E, E_last, n;
	double KAPPASQ = kappa * kappa;
	vector3 r0, v0;
	matrix33 S;
	E_last = M;
	E = M + e * sin(M);
	k = 0;
	while (fabs(E - E_last) > eps && k < NMAX)
	{
		E_last = E;
		E = M + e * sin(E_last);
		k++;
	}
	n = sqrt(KAPPASQ) * pow(a, -1.5);
	r0.x = a * (cos(E) - e);
	r0.y = a * sqrt(1 - e * e) * sin(E);
	r0.z = 0.0;
	v0.x = -n * a * a / norm3(r0) * sin(E);
	v0.y = n * a * a / norm3(r0) * sqrt(1.0 - e * e) * cos(E);
	v0.z = 0.0;
	S = new_matrix(\
		cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i), \
		- cos(Omega) * sin(omega) - sin(Omega) * cos(omega) * cos(i), \
		sin(Omega) * sin(i), \
		sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i), \
		- sin(Omega) * sin(omega) + cos(Omega) * cos(omega) * cos(i), \
		- cos(Omega) * sin(i), \
		sin(omega) * sin(i), \
		cos(omega) * sin(i), \
		cos(i));
	*r = linmap(S, r0);
	*v = linmap(S, v0);
}

void orbital_elements(double kappa, vector3 r, vector3 v, double* a, double* e, double* i, double* Omega, double* omega, double* M)
{
	double KAPPASQ = kappa * kappa;
	vector3 n;
	double c, p, rr, u, theta, dr, E, beta;
	rr = norm3(r);
	*a = 1.0 / (2.0 / rr - dot(v, v) / KAPPASQ);
	n = cross(r, v);
	c = norm3(n);
	n = cdot3(1.0 / c, n);
	*Omega = to_interval(atan2(n.x, -n.y), 0.0, 2.0 * M_PI);
	*i = atan2(sqrt(n.x * n.x + n.y * n.y), n.z);
	p = c * c / KAPPASQ;
	*e = sqrt(1.0 - p / *a);
	u = atan2(r.z / sin(*i), r.x * cos(*Omega) + r.y * sin(*Omega));
	dr = dot(r, v) / rr;
	theta = atan2(dr * sqrt(p) / kappa, p / rr - 1.0);
	*omega = to_interval(u - theta, 0.0, 2.0 * M_PI);
	//beta = *e / (1.0 - (*e) * (*e));
	//E = theta - 2.0 * atan(beta * sin(theta) / (1.0 + beta * cos(theta)));
	E = 2.0 * atan(sqrt((1.0 - *e) / (1.0 + *e)) * tan(0.5 * theta));
	*M = E - *e * sin(E);
}

void calculate_orbital_parameters(int nt, int n_bodies, double kappa, double* coords, double* parameters)
{
	const int N = 6;
	int i, j;
	double* params_cur;
	double* coords_cur;
	double* coords_cur_cur;
	double* params_cur_cur;
	double* a_ptr;
	double* e_ptr;
	double* i_ptr;
	double* Omega_ptr;
	double* omega_ptr;
	double* M_ptr;
	vector3 r, v;
	for (i = 0; i <= nt; i++)
	{
		params_cur = parameters + i * n_bodies * N;
		coords_cur = coords + i * n_bodies * N;
		for (j = 0; j < n_bodies; j++)
		{
			coords_cur_cur = coords_cur + j * N;
			params_cur_cur = params_cur + j * N;
			a_ptr = params_cur_cur;
			e_ptr = params_cur_cur + 1;
			i_ptr = params_cur_cur + 2;
			Omega_ptr = params_cur_cur + 3;
			omega_ptr = params_cur_cur + 4;
			M_ptr = params_cur_cur + 5;
			r = new_vector(coords_cur_cur[0], coords_cur_cur[1], coords_cur_cur[2]);
			v = new_vector(coords_cur_cur[3], coords_cur_cur[4], coords_cur_cur[5]);
			orbital_elements(kappa, r, v, a_ptr, e_ptr, i_ptr, Omega_ptr, omega_ptr, M_ptr);
		}
	}
}

double degrees(double radians)
{
    return radians * 180.0 / M_PI;
}

double radians(double degrees)
{
    return degrees * M_PI / 180.0;
}

double to_interval(double value, double left_bound, double right_bound)
{
	double len, periods;
    len = right_bound - left_bound;
    periods = (value - left_bound) / len;
    return value - floor(periods) * len;
}
