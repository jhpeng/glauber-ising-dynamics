#include <gsl/gsl_rng.h>

#ifndef dtsw_method_h
#define dtsw_method_h

void dtsw_setup(int n, int r, int t_max, int type, double beta, double p, gsl_rng* rng);

void dtsw_update(gsl_rng* rng);

void dtsw_measurement_sampling();

void dtsw_measurement_save();

#endif
