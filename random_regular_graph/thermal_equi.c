#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "gsl/gsl_rng.h"

#include "rrg_generator.h"
#include "dtmc_method.h"

static int nskip=10;

int main(int argc, char** argv) {
    int n=atoi(argv[1]);
    int r=atoi(argv[2]);
    double beta=atof(argv[3]);
    int nsample=atoi(argv[4]);
    int nblock=atoi(argv[5]);
    int nthermal=nsample;
    unsigned long int seed=atoi(argv[6]);
    int nspin=r+1;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);

    // generate random regular graph
    int* graph=NULL;
    int size;

    if(random_regular_graph_generator(n,r,2,&graph,&size,rng))
        return 1;

    while(graph_analysis(n,r,graph)!=1) {
        random_regular_graph_generator(n,r,2,&graph,&size,rng);
    }

    int* state = (int*)malloc(sizeof(int)*n);
    int* state_temp = (int*)malloc(sizeof(int)*n);
    int* temp;

    glauber_ising_transition_prob(nspin,beta);
    dtmc_initial_state(n,0,0.5,state,rng);

    // thermalization
    for(int i=0;i<nthermal;i++) {
        dtmc_update(n,nspin,graph,state,state_temp,rng);
        
        // switch the address between state and state_temp
        temp = state;
        state = state_temp;
        state_temp = temp;
    }

    // measurement
    FILE* m_file = fopen("data.txt","a");
    int bsize = nsample/nblock;
    int Mz;
    double mz,mz1,mz2,mz4;
    for(int i_block=0;i_block<nblock;i_block++) {
        mz1=0;
        mz2=0;
        mz4=0;
        for(int j=0;j<bsize;j++) {
            for(int k=0;k<nskip;k++) {
                dtmc_update(n,nspin,graph,state,state_temp,rng);
        
                // switch the address between state and state_temp
                temp = state;
                state = state_temp;
                state_temp = temp;
            }

            Mz=0;
            for(int i=0;i<n;i++) {
                Mz+=state[i];
            }

            mz = ((double)Mz*2/n-1);
            mz1 += fabs(mz);
            mz2 += mz*mz;
            mz4 += mz*mz*mz*mz;
        }

        mz1 = mz1/bsize;
        mz2 = mz2/bsize;
        mz4 = mz4/bsize;

        fprintf(m_file,"%.12e %.12e %.12e \n",mz1,mz2,mz4);
    }

    fclose(m_file);
}
