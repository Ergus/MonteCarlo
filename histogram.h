/**
 * This is the header file for the histogram in the Montecarlo Exercise
 * Jimmy Aguilar Mena
 * 20/02/2018
 */

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

struct neutron;

//! Histogram struct
struct histogram {
	double len;
	double bin;
	int nbins;
	int *data_count;
	double *data_energy;
	int ncounts, nenergies;
};

//! Function to create histogram struct pointer
struct histogram *create_hist(double sigma, double n, int nbins);

//! Function to release the histogram
void free_hist(struct histogram *in);

//! Register energy in a point for the histograms
int register_energy(struct histogram *in, double energy, struct neutron *part);

//! Register the total len in another histogram
int register_len(struct histogram *in, struct neutron *part);

//! save Histogram to a file
void save_hist(struct histogram *in, const char *filename);

#endif
