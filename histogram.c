/**
 * This is the source file for the histogram in the Montecarlo exercise.
 * Jimmy Aguilar Mena
 * 20/02/2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#include "histogram.h"
#include "neutron.h"


//! Function to create histogram struct pointer
struct histogram *create_hist(double sigma, double n, int nbins)
{
	struct histogram *out = malloc(sizeof(struct histogram));
	assert(out);
	out->ncounts = 0;
	out->nenergies = 0;
	out->len = 10*log(1000) / sigma / n;
	out->nbins = nbins;
	out->bin = out->len / nbins;
	out->data_count = calloc(nbins, sizeof(int));
	out->data_energy = calloc(nbins, sizeof(double));
	assert(out->data_count);
	assert(out->data_energy);
	printf("Histogram len = %lf nbins = %d bin = %lf\n",
	       out->len, out->nbins, out->bin);
	return out;
}

//! Function to release the histogram
void free_hist(struct histogram *in)
{
	free(in->data_count);
	free(in->data_energy);
	free(in);
}

//! Register energy in a point for the histograms
int register_energy(struct histogram *in, double energy, struct neutron *part)
{
	int bin = part->pos[0] / in->bin;
	if (in->nbins <= bin || bin < 0)
		return 0;

	in->data_energy[bin] += energy;
	in->nenergies++;

	return 1;
}

//! Register the total len in another histogram
int register_len(struct histogram *in, struct neutron *part)
{
	int bin = part->pos[0] / in->bin;
	if (in->nbins <= bin || bin < 0)
		return 0;

	in->data_count[bin] ++;
	in->ncounts++;

	return 1;
}

//! save Histogram to a file
void save_hist(struct histogram *in, const char *filename)
{
	time_t timer;
    char buffer[26];
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);


	FILE *out = fopen(filename, "w");
	assert(out);
	int i;
	fprintf(out,"# Run: %s\n",buffer);
	fprintf(out, "# i bin len energy\n");
	fprintf(out, "# ncounts %d nenergies %d\n",
	        in->ncounts, in->nenergies);

	for(i = 0; i < in->nbins; ++i){
		fprintf(out, "%d %lf %d %lf\n",
		        i, (double)i * in->bin, in->data_count[i], in->data_energy[i]);
	}
	fclose(out);
}
