#include <gsl/gsl_rng.h>

#ifndef dtmc_method_h
#define dtmc_method_h

void glauber_ising_transition_prob(int nspin, double beta);

void dtmc_update(int nsite, int nspin, int* graph, int* state_p, int* state_n, gsl_rng* rng);

void dtmc_initial_state(int nsite, int type, double p, int* state, gsl_rng* rng);

int dtmc_final_state(int nsite, int type, int nfix, int* state, gsl_rng* rng);

int dtmc_measurement_count();

void dtmc_measurement_setup(int nsite, int t_max);

void dtmc_measurement_sampling(int nsite, int t_max, int t, int* state);

void dtmc_measurement_converge(int nsite, int t_max);

void dtmc_measurement_save(int nsite, int t_max);

void dtmc_free();

#endif
