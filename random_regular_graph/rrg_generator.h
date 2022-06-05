#include "gsl/gsl_rng.h"

#ifndef rrg_generator_h
#define rrg_generator_h

int random_regular_graph_generator(int n, int r, int m, int* graph, gsl_rng* rng);

int graph_analysis(int n, int r, int* graph);

#endif
