#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "gsl/gsl_rng.h"

#include "rrg_generator.h"
#include "dtmc_method.h"
#include "dtsw_method.h"

int main(int argc, char** argv) {
    int n=atoi(argv[1]);
    int r=atoi(argv[2]);
    double beta=atof(argv[3]);
    double p=atof(argv[4]);
    int nfix=atoi(argv[5]);
    int t_max=atoi(argv[6]);
    int nsample=atoi(argv[7]);
    int nblock=atoi(argv[8]);
    int nthermal=nsample;
    unsigned long int seed=atoi(argv[9]);
    int mode=atoi(argv[10]);
    int nspin=r+1;

    if(mode==0 || mode==1) {
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

    for(int i=0;i<n;i++) {
        printf("i=%d : [",i);
        for(int j=0;j<r;j++) {
            printf("%d ",graph[i*r+j]);
        }
        printf("]\n");
    }


    // setup for dtmc 
    int type=2;
    int ntrial=0;
    int* state = (int*)malloc(sizeof(int)*n);
    int* state_temp = (int*)malloc(sizeof(int)*n);
    int* temp;

    glauber_ising_transition_prob(nspin,beta);

    for(int i=0;i<nblock;i++) {
        dtmc_measurement_setup(n,t_max);

        while(dtmc_measurement_count()<(nsample/nblock)) {

            dtmc_initial_state(n,type,p,state,rng);
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
            if(dtmc_final_state(n,type,nfix,state,rng)) {
                dtmc_measurement_converge(n,t_max);

                printf("samples %d/%d | nblock %d/%d",dtmc_measurement_count(),(nsample/nblock),i,nblock);
                printf("\r");
                fflush(stdout);
            }
            ntrial++;
        }

        dtmc_measurement_save(n,t_max);
    }
    printf("samples %d/%d | nblock %d/%d\n",(nsample/nblock),(nsample/nblock),nblock,nblock);

    // save information about dtmc
    FILE* dtmc_file=fopen("dtmc_info.txt","w");

    fprintf(dtmc_file,"#dtmc_info.txt");
    fprintf(dtmc_file,"n      : %d\n",n);
    fprintf(dtmc_file,"r      : %d\n",r);
    fprintf(dtmc_file,"beta   : %.6f\n",p);
    fprintf(dtmc_file,"p      : %.6f\n",p);
    fprintf(dtmc_file,"nfix   : %d\n",nfix);
    fprintf(dtmc_file,"t_max  : %d\n",t_max);
    fprintf(dtmc_file,"nsample: %d\n",nsample);
    fprintf(dtmc_file,"ntrial : %d\n",ntrial);
    fclose(dtmc_file);

    }
    
    if(mode==0 || mode==2) {
    // dtsw_method
    // type=0 : fixed nfix final sites
    // type=1 : fixed all initial and final sites
    int type=1;
    printf("===================================\n");
    printf("starting dtsw method...\n");
    gsl_rng* rng2 = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng2, seed);

    printf("allocating memory...\n");
    dtsw_setup(n,r,t_max,type,nfix,beta,p,rng2);
    
    printf("thermalization...\n");
    for(int i=0;i<nthermal;i++) {
        dtsw_update(rng2);
        printf("thernalization %d/%d",i,nthermal);
        printf("\r");
        fflush(stdout);
    }
    printf("thernalization %d/%d\n",nthermal,nthermal);

    printf("starting the estimator...\n");
    for(int i=0;i<nblock;i++) {
        for(int j=0;j<(nsample/nblock);j++) {
            for(int k=0;k<10;k++) dtsw_update(rng2);
            dtsw_measurement_sampling();
        }
        dtsw_measurement_save();
    }

    }

    return 0;
}
