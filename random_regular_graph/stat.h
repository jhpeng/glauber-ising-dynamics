#ifndef stat_h
#define stat_h


void stat_setup(int ndata, int nbin);

void stat_append(double* data);

void stat_print_status();

void stat_free();

#endif
