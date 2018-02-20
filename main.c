/**
 * This is the source file for the main in the Montecarlo Exercise
 * Jimmy Aguilar Mena
 * 20/02/2018
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "histogram.h"
#include "neutron.h"
#include "functions.h"

 //! All the magnitudes in meters
double cut_energy = 1E-3;

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
	double atom_mass = get_mass(A);
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
			.m = get_mass(1),
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
