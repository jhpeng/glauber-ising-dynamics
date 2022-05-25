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
    } else if(type==1) {
        int total_m=2;
        while(total_m>1 || total_m<-1) {
            total_m=0;
            for(int i=0;i<nsite;i++) {
                state[i]=0;
                if(gsl_rng_uniform_pos(rng)<p) 
                    state[i]=1;

                total_m += state[i]*2-1;
            }
        }
    }
}

int dtmc_final_state(int nsite, int type, int nfix, int* state, gsl_rng* rng) {
    int check=1;
    if(type==0) {
        check=1;
        for(int j=0;j<nfix;j++) {
            if(state[j]==0) check=0;
        }
    }

    return check;
}

static int dtmc_count=0;
static double* dtmc_mz_total=NULL;
static double* dtmc_mz_abs_total=NULL;
static double* dtmc_mz_1_one_total=NULL;
static double* dtmc_mz_1_two_total=NULL;
static double* dtmc_mz_hold=NULL;
static double* dtmc_mz_1_hold=NULL;

int dtmc_measurement_count() {
    return dtmc_count;
}

void dtmc_measurement_setup(int nsite, int t_max) {
    if(dtmc_mz_total!=NULL) free(dtmc_mz_total);
    if(dtmc_mz_abs_total!=NULL) free(dtmc_mz_abs_total);
    if(dtmc_mz_1_one_total!=NULL) free(dtmc_mz_1_one_total);
    if(dtmc_mz_1_two_total!=NULL) free(dtmc_mz_1_two_total);
    if(dtmc_mz_hold!=NULL) free(dtmc_mz_hold);
    if(dtmc_mz_1_hold!=NULL) free(dtmc_mz_1_hold);

    dtmc_mz_total = (double*)malloc(sizeof(double)*(t_max+1));
    dtmc_mz_abs_total = (double*)malloc(sizeof(double)*(t_max+1));
    dtmc_mz_1_one_total = (double*)malloc(sizeof(double)*(t_max+1));
    dtmc_mz_1_two_total = (double*)malloc(sizeof(double)*(t_max+1)*(t_max+1));
    dtmc_mz_hold = (double*)malloc(sizeof(double)*(t_max+1));
    dtmc_mz_1_hold = (double*)malloc(sizeof(double)*(t_max+1));

    for(int t=0;t<(t_max+1);t++) {
        dtmc_mz_total[t]=0;
        dtmc_mz_abs_total[t]=0;
        dtmc_mz_hold[t]=0;
        dtmc_mz_1_hold[t]=0;
        dtmc_mz_1_one_total[t]=0;

        for(int s=0;s<(t_max+1);s++) {
            dtmc_mz_1_two_total[t*(t_max+1)+s]=0;
        }
    }

    dtmc_count=0;
}

void dtmc_measurement_sampling(int nsite, int t_max, int t, int* state) {
    dtmc_mz_hold[t]=0;

    for(int i=0;i<nsite;i++) {
        dtmc_mz_hold[t] += state[i];
    }

    dtmc_mz_hold[t] = 2*dtmc_mz_hold[t]/nsite-1;
    dtmc_mz_1_hold[t] = state[0]*2-1;
}

void dtmc_measurement_converge(int nsite, int t_max) {
    
    for(int t=0;t<(t_max+1);t++) {
        dtmc_mz_total[t] += dtmc_mz_hold[t];
        dtmc_mz_abs_total[t] += fabs(dtmc_mz_hold[t]);

        dtmc_mz_1_one_total[t] += dtmc_mz_1_hold[t];
        for(int s=0;s<(t_max+1);s++) {
            double mz_t = dtmc_mz_1_hold[t];
            double mz_s = dtmc_mz_1_hold[s];
            dtmc_mz_1_two_total[t*(t_max+1)+s] += mz_t*mz_s;
        }
    }

    dtmc_count++;
}

void dtmc_measurement_save(int nsite, int t_max) {
    FILE* mz_file = fopen("dtmc_mz.txt","a");
    FILE* mz_abs_file = fopen("dtmc_mz_abs.txt","a");
    FILE* mz_1_one_file = fopen("dtmc_mz_1_one.txt","a");
    FILE* mz_1_two_file = fopen("dtmc_mz_1_two.txt","a");

    for(int t=0;t<(t_max+1);t++) {
        fprintf(mz_file,"%.16e ", dtmc_mz_total[t]/dtmc_count);
        fprintf(mz_abs_file,"%.16e ", dtmc_mz_abs_total[t]/dtmc_count);
        fprintf(mz_1_one_file,"%.16e ", dtmc_mz_1_one_total[t]/dtmc_count);
        for(int s=0;s<(t_max+1);s++) {
            fprintf(mz_1_two_file,"%.16e ", dtmc_mz_1_two_total[t*(t_max+1)+s]/dtmc_count);
        }
    }
    fprintf(mz_file,"\n");
    fprintf(mz_abs_file,"\n");
    fprintf(mz_1_one_file,"\n");
    fprintf(mz_1_two_file,"\n");

    fclose(mz_file);
    fclose(mz_abs_file);
    fclose(mz_1_one_file);
    fclose(mz_1_two_file);
}
