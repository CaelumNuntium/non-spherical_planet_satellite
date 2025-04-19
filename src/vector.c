#include <stdio.h>
#include <math.h>
#include "vector.h"

vector3 new_vector(double x, double y, double z)
{
	vector3 res;
	res.x = x;
	res.y = y;
	res.z = z;
	return res;
}

double norm3(vector3 v)
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double dot(vector3 a, vector3 b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

vector3 cdot3(double a, vector3 v)
{
	vector3 res;
	res.x = v.x * a;
	res.y = v.y * a;
	res.z = v.z * a;
	return res;
}

vector3 sum(vector3 a, vector3 b)
{
	vector3 res;
	res.x = a.x + b.x;
	res.y = a.y + b.y;
	res.z = a.z + b.z;
	return res;
}

vector3 negative(vector3 v)
{
	vector3 res;
	res.x = -v.x;
	res.y = -v.y;
	res.z = -v.z;
	return res;
}

vector3 cross(vector3 a, vector3 b)
{
	vector3 res;
	res.x = a.y * b.z - a.z * b.y;
	res.y = a.z * b.x - a.x * b.z;
	res.z = a.x * b.y - a.y * b.x;
	return res;
}

void print_vector(vector3 v)
{
	printf("(%lf, %lf, %lf)\n", v.x, v.y, v.z);
}

matrix33 new_matrix(double m00, double m01, double m02, double m10, double m11, double m12, double m20, double m21, double m22)
{
	matrix33 res;
	res.xx = m00;
	res.xy = m01;
	res.xz = m02;
	res.yx = m10;
	res.yy = m11;
	res.yz = m12;
	res.zx = m20;
	res.zy = m21;
	res.zz = m22;
	return res;
}

double get_element(matrix33 m, int row, int col)
{
	if (row > 2 || col > 2 || row < 0 || col < 0)
	{
		fprintf(stderr, "Error: Invalid index (%d, %d) in matrix 3 x 3\n", row, col);
		return NAN;
	}
	switch (row * 3 + col)
	{
	case 0:
		return m.xx;
	case 1:
		return m.xy;
	case 2:
		return m.xz;
	case 3:
		return m.yx;
	case 4:
		return m.yy;
	case 5:
		return m.yz;
	case 6:
		return m.zx;
	case 7:
		return m.zy;
	case 8:
		return m.zz;
	default:
		return NAN;
	}
}

vector3 get_row(matrix33 m, int i)
{
	vector3 res;
	if (i < 0 || i > 2)
	{
		fprintf(stderr, "Error: Invalid index (%d, :) in matrix 3 x 3", i);
		res.x = NAN;
		res.y = NAN;
		res.z = NAN;
	}
	res.x = get_element(m, i, 0);
	res.y = get_element(m, i, 1);
	res.z = get_element(m, i, 2);
	return res;
}

vector3 get_column(matrix33 m, int i)
{
	vector3 res;
	if (i < 0 || i > 2)
	{
		fprintf(stderr, "Error: Invalid index (:, %d) in matrix 3 x 3", i);
		res.x = NAN;
		res.y = NAN;
		res.z = NAN;
	}
	res.x = get_element(m, 0, i);
	res.y = get_element(m, 1, i);
	res.z = get_element(m, 2, i);
	return res;
}

vector3 linmap(matrix33 m, vector3 v)
{
	vector3 res;
	res.x = dot(get_row(m, 0), v);
	res.y = dot(get_row(m, 1), v);
	res.z = dot(get_row(m, 2), v);
	return res;
}

void vector_to_array(vector3 v, double* arr)
{
	arr[0] = v.x;
	arr[1] = v.y;
	arr[2] = v.z;
}

void matrix_to_array(matrix33 m, double* arr)
{
	arr[0] = m.xx;
	arr[1] = m.xy;
	arr[2] = m.xz;
	arr[3] = m.yx;
	arr[4] = m.yy;
	arr[5] = m.yz;
	arr[6] = m.zx;
	arr[7] = m.zy;
	arr[8] = m.zz;
}
