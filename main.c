/**
 * This is the source file for the main in the Montecarlo Exercise
 * Jimmy Aguilar Mena
 * 20/02/2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "histogram.h"
#include "neutron.h"

 //! All the magnitudes in meters
double r_0 = 1.25E-15;
double pi = 3.14;
double unit_mass = 1.6750E-27;
double cut_energy = 1E-3;

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

//! \param[in] y Uniform (0;1] random number
//! \param[in] sigma Cross Section
//! \param[in] n Number of nucleus/volume units.
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

int main(int argc, char** argv)
{
	if (argc < 5){
		printf("Usage: %s A_material Neutron_Energy density events\n", argv[0]);
		return -1;
	}

	int i;
	int A = atoi(argv[1]);
	double E_0 = atof(argv[2]);
	double density = atof(argv[3]);
	int events = atoi(argv[4]);

	double sigma = get_sigma(A);
	double n = get_n(density, A);
	double atom_mass = A * unit_mass;
	int nprint = events / 10;

	printf("Start Run\n");
	printf("A = %d E_0 = %lf MeV density = %lf events %d\n",
	       A, E_0, density, events);

	struct histogram *hist = create_hist(sigma, n, 1000);

	for (i = 0; i < events; ++i) {
		struct neutron in = {
			.v = {1., 0., 0.},
			.ang = {0.,0.},
			.pos = {0., 0., 0.},
			.m = unit_mass,
			.energy = E_0
		};

		in.energy = E_0;
		in.v[0] = 1.0;

		while (in.energy > cut_energy) {
			// Here starts;
			double len = l(rand_uniform(), sigma, n);

			double dx = len * in.v[0];
			double dy = len * in.v[1];
			double dz = len * in.v[2];

			// relative component angles
			double ang0 = get_ang1();
			double ang1 = get_ang2(ang1);

			// relative total angle
			double alpha = get_alpha(ang0, ang1);
			double Ef = get_Ef(alpha, &in, atom_mass);

			assert(Ef < in.energy || "Energy conservation check");

			double DE = in.energy - Ef;  // Deposited energy

			// Set values in the struct
			in.energy = Ef;

			// Absolute final angle
			in.ang[0] += ang0;
			in.ang[1] += ang1;

			in.v[0] = cos(ang0)*cos(ang1);
			in.v[1] = cos(ang0)*sin(ang1);
			in.v[2] = sin(ang0);

			in.pos[0] += dx;
			in.pos[1] += dy;
			in.pos[2] += dz;

			int registered = register_energy(hist, DE, &in);
			if (!registered)
				printf("Error: %d Hit: %lf %lf %lf angle %lf E_0 = %lf Ef = %lf\n",
				       i, in.pos[0], in.pos[1], in.pos[2], alpha, in.energy, Ef);
		}
		register_energy(hist, in.energy, &in);
		register_len(hist, &in);

		if(! (i % nprint))
			printf("%d: Pos: %lf %lf %lf angle %lf Ef = %lf\n",
			       i, in.pos[0], in.pos[1], in.pos[2], in.energy);

	}

	save_hist(hist, "output.out");
	printf("Saved histogram\n");
	free_hist(hist);
	return 0;
}
