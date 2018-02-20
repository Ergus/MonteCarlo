#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "neutron.h"

//! These are constants, you could improve this with macros.
double r_0 = 1.25E-15;
double pi = 3.14;
double unit_mass = 1.6750E-27;

//! Returns the atom mass
double get_mass(int A)
{
	return A * unit_mass;
}

//! This generates a random number in (0;1]
double rand_uniform()
{
	return (double)rand() / (double)((unsigned)RAND_MAX + 1);
}

//! This generates a random number in [-1;1]
double rand_simetric()
{
	return 2 * (double)rand() / (double)((unsigned) RAND_MAX) - 1;
}

//! This generate a random with the step distribution.
double l(double y, double sigma, double n)
{
	return log(1 / (1 - y)) / n / sigma;
}

//! Calculated the cross section
double get_sigma(int A)
{
	double R =  pow((double)A, 1/3) + 1;
	return pi * pow(r_0, 2) * pow(R, 2);
}

//! Calculate relative zenital angle
double get_ang1()
{
	double x = rand_simetric();
	return asin(x);
}

//! calculate relative asimutal angle
double get_ang2(double ang1)
{
	double x = sin(ang1);
	double y = rand_simetric() * sqrt( 1 - pow(x,2) );
	return asin(y);
}

//! Calculates the total impact angle
double get_alpha(double ang1, double ang2)
{
	double c1 = cos(ang1);
	double s2 = sin(ang2);
	double c2 = cos(ang2);

	double denom = pow(sin(ang2),2) - pow(sin(ang2),2) / pow(c1, 2) + 1 / pow(c1,2);
	double cos_alpha = sqrt( pow(c2,2) / denom );

	return acos(cos_alpha);
}

//! Calculates the final energy
double get_Ef(double alpha, struct neutron *in, double atom_mass)
{
	double E = in->energy;
	double a = in->m / atom_mass ;
	double frac = (cos(alpha) - a) / (cos(alpha) + a);
	return E * pow(frac, 2);
}

double get_n(int density, double A)
{
	double mass_atom = A * unit_mass;
	return density / mass_atom;
}

//! Returns speed from energy
double to_speed(double energy)
{
	return sqrt(2 * energy / unit_mass);
}
