#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define TEST_MODE

int* sigma;
int* sigma_temp;

int Nfix = 0;
double* local_prob;
double IAM=0;
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
    } else if(type==3) {
        for(int i=0;i<nsite;i++) {
            sigma[i] = 0;
            if(gsl_rng_uniform_pos(rng)<IAM) {
                sigma[i] = 1;
            }
        }
    } else if(type==4) {
        for(int i=0;i<nsite;i++) {
            sigma[i] = 0;
            if(gsl_rng_uniform_pos(rng)<IAM || i<Nfix) {
                sigma[i] = 1;
            }
        }
    }
}

int check_final_state(int nsite) {
    int check=1;
    for(int i=0;i<Nfix;i++) {
        if(sigma[i]==0)
            check=0;
    }

    return check;
}

int* RandomGraphConf;
void generate_random_graph(int nsite, unsigned long int seed) {
    if(RandomGraphConf==NULL) {
        RandomGraphConf = (int*)malloc(sizeof(int)*nsite*2);
    } else {
        free(RandomGraphConf);
        RandomGraphConf = (int*)malloc(sizeof(int)*nsite*2);
    }

    gsl_rng* rng_temp = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng_temp,seed);

    for(int x=0;x<nsite;x++) {
        int x1 = (int)(gsl_rng_uniform_pos(rng_temp)*nsite);
        int x2 = (int)(gsl_rng_uniform_pos(rng_temp)*nsite);

        while(x1==x2 || x1==x || x2==x) {
            x1 = (int)(gsl_rng_uniform_pos(rng_temp)*nsite);
            x2 = (int)(gsl_rng_uniform_pos(rng_temp)*nsite);
        }
        
        RandomGraphConf[x*2+0] = x1;
        RandomGraphConf[x*2+1] = x2;
    }
    gsl_rng_free(rng_temp);
}

double WeightOfConf;
void update(int nsite) {
    int s;
    for(int i=0;i<nsite;i++) {
        s = sigma[i]+sigma[RandomGraphConf[i*2+0]]+sigma[RandomGraphConf[i*2+1]];
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

#ifdef TEST_MODE
unsigned long int ITERATION_TOTAL=0;
#endif

int main(int argc, char** argv) {

    int nsite = atoi(argv[1]);
    int t = atoi(argv[2]);
    //double mz = 0;
    double mz = atof(argv[3]);
    IAM = (mz+1)*0.5;
    Nfix = atoi(argv[4]);
    double beta = atof(argv[5]);
    int nblock = atoi(argv[6]);
    int nsample = atoi(argv[7]);
    unsigned long int seed = atoi(argv[8]);

    generate_random_graph(nsite,seed);

    sigma      = (int*)malloc(sizeof(int)*nsite);
    sigma_temp = (int*)malloc(sizeof(int)*nsite);

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    set_local_prob(beta);

    double* magz = (double*)malloc(sizeof(double)*t);
    double* magz_temp = (double*)malloc(sizeof(double)*t);

    char magz_filename[128];
    char weight_filename[128];
    //sprintf(magz_filename,"magz_l_%d_mz_%.6f_beta_%.6f_t_%d.txt",nsite,mz,beta,t);
    //sprintf(weight_filename,"weight_l_%d_mz_%.6f_beta_%.6f_t_%d.txt",nsite,mz,beta,t);
    sprintf(magz_filename,"magz_l_%d_nfix_%d_beta_%.6f_t_%d_seed_%ld_.txt",nsite,Nfix,beta,t,seed);
    sprintf(weight_filename,"weight_l_%d_nfix_%d_beta_%.6f_t_%d_seed_%ld_.txt",nsite,Nfix,beta,t,seed);
    FILE* magz_file = fopen(magz_filename,"w");
    FILE* weight_file = fopen(weight_filename,"w");

    for(int k=0;k<nblock;k++) {
        for(int i=0;i<t;i++) magz[i]=0;
        for(int j=0;j<nsample;j++) {
            int check=0;
            while(check==0){
                WeightOfConf = 0;
                initial_state(4,nsite);
                for(int i=0;i<t;i++) {
                    //show_state(nsite);
                    update(nsite);
                    magz_temp[i]=0;
                    for(int i_site=0;i_site<nsite;i_site++) {
                        magz_temp[i] += sigma[i_site];
                    }
                }
                check = check_final_state(nsite);

#ifdef TEST_MODE
                ITERATION_TOTAL++;
#endif
            }

            for(int i=0;i<t;i++) 
                magz[i]+=magz_temp[i];

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

#ifdef TEST_MODE
    printf("total iteration : %ld\n",ITERATION_TOTAL);
    FILE* fite = fopen("total_iteration.txt","a");
    fprintf(fite, "%d %d %lf %d %lf %d %d %ld %ld\n",nsite,t,mz,Nfix,beta,nblock,nsample,seed,ITERATION_TOTAL);
#endif

    return 0;
}
