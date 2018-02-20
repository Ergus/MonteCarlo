#ifndef FUNCTIONS_H
#define FUNCTIONS_H

//! This generates a random number in (0;1]
double rand_uniform();

//! This generates a random number in [-1;1]
double rand_simetric();

//! Returns the atom mass
double get_mass(int A);

//! This generate a random with the step distribution.

//! \param[in] y Uniform (0;1] random number
//! \param[in] sigma Cross Section
//! \param[in] n Number of nucleus/volume units.
double l(double y, double sigma, double n);

//! Calculated the cross section
double get_sigma(int A);

//! Calculate relative zenital angle
double get_ang1();

//! calculate relative asimutal angle
double get_ang2(double ang1);

//! Calculates the total impact angle
double get_alpha(double ang1, double ang2);

//! Calculates the final energy
double get_Ef(double alpha, struct neutron *in, double atom_mass);

//! Returns the number of nucleous / distance unit
double get_n(int density, double A);

//! Returns speed from energy
double to_speed(double energy);

#endif
