#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int* sigma;
int* sigma_temp;

double* local_prob;

void update(int** sigma, int** placeholder, double beta, int nsite) {
    if(local_prob==NULL) {
        local_prob = (double*)malloc(sizeof(double)*4);
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

    int s,ss;
    for(int i=0;i<nsite;i++) {
        s = (*sigma)[i]+(*sigma)[(i+1)%nsite]+(*sigma)[(i-1+nsite)%nsite];
    }
}

int main() {

    int N = 32;
    double Beta = 0.5;
    sigma      = (int*)malloc(sizeof(int)*N);
    sigma_temp = (int*)malloc(sizeof(int)*N);

    update(&sigma,&sigma_temp,Beta,N);

    return 0;
}
