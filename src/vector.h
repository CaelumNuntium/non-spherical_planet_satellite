#pragma once

typedef struct {
	double x, y, z;
} vector3;

typedef struct {
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;
} matrix33;

double norm3(vector3);
double dot(vector3, vector3);
vector3 cdot3(double, vector3);
vector3 sum(vector3, vector3);
vector3 negative(vector3);
vector3 new_vector(double, double, double);
vector3 cross(vector3, vector3);
void print_vector(vector3);
matrix33 new_matrix(double, double, double, double, double, double, double, double, double);
vector3 linmap(matrix33, vector3);
vector3 get_row(matrix33, int);
vector3 get_column(matrix33, int);
double get_element(matrix33, int, int);
void vector_to_array(vector3, double*);
void matrix_to_array(matrix33, double*);
