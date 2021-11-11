#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int* sigma;
int* sigma_temp;

double* local_prob;
gsl_rng* rng;


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
    } else if(type==2) {
        for(int i=0;i<nsite;i++) {
            sigma[i] = 1;
            if(i%4==0) sigma[i] = 0;
        }
    }
}

double WeightOfConf;
void update(int nsite) {
    int s;
    for(int i=0;i<nsite;i++) {
        s = sigma[i]+sigma[(i+1)%nsite]+sigma[(i-1+nsite)%nsite];
        if(gsl_rng_uniform_pos(rng)<local_prob[s]) {
            sigma_temp[i]=1;
            WeightOfConf += log(local_prob[s]);
        } else { 
            sigma_temp[i]=0;
            WeightOfConf += log(1-local_prob[s]);
        }
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
    int nblock = atoi(argv[4]);
    int nsample = atoi(argv[5]);
    unsigned long int seed = atoi(argv[6]);

    sigma      = (int*)malloc(sizeof(int)*nsite);
    sigma_temp = (int*)malloc(sizeof(int)*nsite);

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    set_local_prob(beta);

    double* magz = (double*)malloc(sizeof(double)*t);

    char magz_filename[128];
    char weight_filename[128];
    sprintf(magz_filename,"magz_l_%d_beta_%.6f_t_%d.txt",nsite,beta,t);
    sprintf(weight_filename,"weight_l_%d_beta_%.6f_t_%d.txt",nsite,beta,t);
    FILE* magz_file = fopen(magz_filename,"w");
    FILE* weight_file = fopen(weight_filename,"w");

    for(int k=0;k<nblock;k++) {
        for(int i=0;i<t;i++) magz[i]=0;
        for(int j=0;j<nsample;j++) {
            WeightOfConf = 0;
            initial_state(2,nsite);
            for(int i=0;i<t;i++) {
                //show_state(nsite);
                update(nsite);
                for(int i_site=0;i_site<nsite;i_site++) {
                    magz[i]+=sigma[i_site];
                }
            }
            if(j%10==0)
                fprintf(weight_file,"%16e ",WeightOfConf);
        }
        fprintf(weight_file,"\n");

        for(int i=0;i<t;i++) {
            magz[i] = magz[i]/nsample/nsite;
            magz[i] = magz[i]*2-1;

            fprintf(magz_file,"%.16e ",magz[i]);
        }
        fprintf(magz_file,"\n");

    }

    fclose(magz_file);
    fclose(weight_file);

    return 0;
}
