#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "gsl/gsl_rng.h"

#include "rrg_generator.h"
#include "dtmc_method.h"
#include "dtsw_method.h"

int main() {
    int n=20;
    int r=3;
    double beta=0.25;
    double p=0.6;
    int t_max=20;
    int nsample=1000000;
    int nblock=1;
    int nthermal=100000;
    unsigned long int seed=84893;
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
/*
    int max_ncluster=0;
    int* ncluster_collect=NULL;
    double total_sampling=0;
    while(1) {
        random_regular_graph_generator(n,r,2,&graph,&size,rng);
        int ncluster = graph_analysis(n,r,graph);
        if(ncluster>max_ncluster) {
            int* temp = (int*)malloc(sizeof(int)*ncluster);
            for(int i=0;i<ncluster;i++) temp[i]=0;
            for(int i=0;i<max_ncluster;i++)
                temp[i]=ncluster_collect[i];

            if(ncluster_collect!=NULL) free(ncluster_collect);
            ncluster_collect=temp;
            max_ncluster = ncluster;
        }
        ncluster_collect[ncluster-1]++;

        total_sampling++;
        for(int i=0;i<max_ncluster;i++)
            printf("total sampling on %d cluster : %d (%.6f)\n",i+1,ncluster_collect[i],ncluster_collect[i]/total_sampling);
    }
*/
    for(int i=0;i<n;i++) {
        printf("i=%d : [",i);
        for(int j=0;j<r;j++) {
            printf("%d ",graph[i*r+j]);
        }
        printf("]\n");
    }


    // setup for dtmc 
    int* state = (int*)malloc(sizeof(int)*n);
    int* state_temp = (int*)malloc(sizeof(int)*n);
    int* temp;

    glauber_ising_transition_prob(nspin,beta);

    for(int i=0;i<nblock;i++) {
        dtmc_measurement_setup(n,t_max);

        while(dtmc_measurement_count()<(nsample/nblock)) {

            dtmc_initial_state(n,0,p,state,rng);
            dtmc_measurement_sampling(n,t_max,0,state);

            for(int t=0;t<t_max;t++) {
                //for(int i=0;i<n;i++) printf("%d ",state[i]);
                //printf("\n");
                dtmc_update(n,nspin,graph,state,state_temp,rng);

                // switch the address between state and state_temp
                temp = state;
                state = state_temp;
                state_temp = temp;

                dtmc_measurement_sampling(n,t_max,t+1,state);
            }
    
            // condition of final state
            dtmc_measurement_converge(n,t_max);
        }

        dtmc_measurement_save(n,t_max);
    }



    // dtsw_method
    //gsl_rng_free(rng);
    gsl_rng* rng2 = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng2, seed);

    dtsw_setup(n,r,t_max,0,beta,p,rng2);
    
    for(int i=0;i<nthermal;i++) {
        dtsw_update(rng);
    }
    for(int i=0;i<nblock;i++) {
        for(int j=0;j<(nsample/nblock);j++) {
            for(int k=0;k<10;k++) dtsw_update(rng);
            dtsw_measurement_sampling();
        }
        dtsw_measurement_save();
    }


    return 0;
}
