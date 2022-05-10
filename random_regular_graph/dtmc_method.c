#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

static double* trans_probs=NULL;

void glauber_ising_transition_prob(int nspin, double beta) {
    double* trans_probs_temp = (double*)malloc(sizeof(double)*(nspin+1));

    double w0,w1,z;
    int total_spin;
    for(int s=0;s<(nspin+1);s++) {
        total_spin = s*2-nspin;
        w1 = exp(beta*total_spin);
        w0 = exp(-beta*total_spin);
        z = w1+w0;
        trans_probs_temp[s] = w1/z;
    }

    if(trans_probs!=NULL)
        free(trans_probs);

    trans_probs = trans_probs_temp;
}

void dtmc_update(int nsite, int nspin, int* graph, int* state_p, int* state_n, gsl_rng* rng) {
    int r=nspin-1; // number of edge for each vertex
    int s;
    for(int i=0;i<nsite;i++) {
        s = state_p[i];
        for(int j=0;j<r;j++)
            s += state_p[graph[i*r+j]];

        if(gsl_rng_uniform_pos(rng)<trans_probs[s]) {
            state_n[i]=1;
        } else {
            state_n[i]=0;
        }
    }
}

void dtmc_initial_state(int nsite, int type, double p, int* state, gsl_rng* rng) {
    if(type==0) {
        for(int i=0;i<nsite;i++) {
            if(gsl_rng_uniform_pos(rng)<p) {
                state[i]=1;
            } else {
                state[i]=0;
            }
        }
    }
}

static int dtmc_count=0;
static double* dtmc_mz_total=NULL;
static double* dtmc_mz_hold=NULL;

int dtmc_measurement_count() {
    return dtmc_count;
}

void dtmc_measurement_setup(int nsite, int t_max) {
    if(dtmc_mz_total!=NULL) free(dtmc_mz_total);
    if(dtmc_mz_hold!=NULL) free(dtmc_mz_hold);

    dtmc_mz_total = (double*)malloc(sizeof(double)*(t_max+1));
    dtmc_mz_hold = (double*)malloc(sizeof(double)*(t_max+1));

    for(int i=0;i<(t_max+1);i++) {
        dtmc_mz_total[i]=0;
        dtmc_mz_hold[i]=0;
    }

    dtmc_count=0;
}

void dtmc_measurement_sampling(int nsite, int t_max, int t, int* state) {
    dtmc_mz_hold[t]=0;

    for(int i=0;i<nsite;i++) {
        dtmc_mz_hold[t] += state[i];
    }

    dtmc_mz_hold[t] = 2*dtmc_mz_hold[t]/nsite-1;
}

void dtmc_measurement_converge(int nsite, int t_max) {
    
    for(int t=0;t<(t_max+1);t++) {
        dtmc_mz_total[t] += dtmc_mz_hold[t];
    }

    dtmc_count++;
}

void dtmc_measurement_save(int nsite, int t_max) {
    FILE* mz_file = fopen("dtmc_mz.txt","a");

    for(int t=0;t<(t_max+1);t++) {
        fprintf(mz_file,"%.16e ", dtmc_mz_total[t]/dtmc_count);
    }
    fprintf(mz_file,"\n");

    fclose(mz_file);
}
