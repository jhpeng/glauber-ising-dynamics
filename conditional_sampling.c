#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int* sigma;
int* sigma_temp;

double* local_prob;
gsl_rng* rng;

double* magz;

void set_local_prob(double beta){
    if(local_prob==NULL) {
        local_prob = (double*)malloc(sizeof(double)*4);
    }
    double p1,p0,z;
    int ss;
    for(int s=0;s<4;s++){
        ss = s*2-3;
        p1 = exp(beta*ss);
        p0 = exp(-beta*ss);
        z=p1+p0;
        p1 = p1/z;
        local_prob[s] = p1;
    }
}

void initial_state(int type, int nsite) {
    if(type==0){
        for(int i=0;i<nsite;i++) {
            sigma[i] = i%2;
        }
    } else if(type==1) {
        for(int i=0;i<nsite;i++) {
            sigma[i] = 1;
        }
    }
}

void update(int nsite) {
    int s;
    for(int i=0;i<nsite;i++) {
        s = sigma[i]+sigma[(i+1)%nsite]+sigma[(i-1+nsite)%nsite];
        if(gsl_rng_uniform_pos(rng)<local_prob[s]) sigma_temp[i]=1;
        else sigma_temp[i]=0;
    }

    int* p = sigma;
    sigma = sigma_temp;
    sigma_temp = p;
}

void show_state(int nsite) {
    for(int i=0;i<nsite;i++) {
        printf("%d ",sigma[i]);
    }
    printf("\n");
}

int main(int argc, char** argv) {

    int nsite = atoi(argv[1]);
    int t = atoi(argv[2]);
    double beta = atof(argv[3]);
    int nsample = atoi(argv[4]);
    unsigned long int seed = atoi(argv[5]);

    sigma      = (int*)malloc(sizeof(int)*nsite);
    sigma_temp = (int*)malloc(sizeof(int)*nsite);

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    set_local_prob(beta);

    magz = (double*)malloc(sizeof(double)*t);
    for(int i=0;i<t;i++) magz[i]=0;

    for(int j=0;j<nsample;j++) {
        initial_state(1,nsite);
        for(int i=0;i<t;i++) {
            //show_state(nsite);
            update(nsite);
            for(int i_site=0;i_site<nsite;i_site++) {
                magz[i]+=sigma[i_site];
            }
        }
    }

    for(int i=0;i<t;i++) {
        magz[i] = magz[i]/nsample/nsite;
        magz[i] = magz[i]*2-1;

        printf("t=%d %.10e\n",i+1,magz[i]);
    }

    return 0;
}
