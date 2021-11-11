#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "union_find.h"

int* Sigma;
int* Bond2index;
int Nsite,Nbond;

double Prob[16];
double LocalWeight[16];

int* Ptree;
int* Weight;
int* Flist;

gsl_rng* rng;

int bits2int(int *b) {
    int i=0;
    i+=b[0];
    i+=(b[1])<<1;
    i+=(b[2])<<2;
    i+=(b[3])<<3;

    return i;
}

void int2bits(int i, int* b) {
    b[0] = i&1;
    i = i>>1;
    b[1] = i&1;
    i = i>>1;
    b[2] = i&1;
    i = i>>1;
    b[3] = i&1;
}

double configuration2weight(int i, double beta) {
    int b[4],s;

    int2bits(i,b);

    s = b[1]+b[2]+b[3];
    s = s*2-3;
    double w = exp((b[0]*2-1)*s*beta);

    return w/(w+1.0/w);
}

void initial_state(int L, int T, int type) {
    if(type==0){
        for(int i=0;i<L;i++) Sigma[i]=i%2;
        for(int i=L;i<L*T;i++) {
            Sigma[i]=1;
            if(gsl_rng_uniform_pos(rng)<0.5){
                Sigma[i]=0;
            }
        }
    } else if(type==1) {
        for(int i=0;i<L;i++) Sigma[i]=1;
        for(int i=L;i<L*T;i++) {
            Sigma[i]=1;
            if(gsl_rng_uniform_pos(rng)<0.5){
                Sigma[i]=0;
            }
        }
    } else if(type==2) {
        for(int i=0;i<L;i++) {
            Sigma[i]=1;
            if(i%4==0) Sigma[i]=0;
        }
        for(int i=L;i<L*T;i++) {
            Sigma[i]=1;
            if(gsl_rng_uniform_pos(rng)<0.5){
                Sigma[i]=0;
            }
        }
    }
}

void initial_tree() {
    for(int i=0;i<Nsite;i++) {
        Ptree[i]  = i;
        Weight[i] = 1;
    }
}

void boundary_condition(int L, int T, int type) {
    if(type==0) {
        // fix the initial state
        // free the final state
        for(int i=0;i<L;i++) {
            Weight[i] = -1;
        }
    } else if(type==1) {
        // fix the initial state
        // free the final state
        for(int i=0;i<L;i++) {
            Weight[i] = -1;
        }
        
        for(int i=(T-1)*L; i<Nsite; i++) {
            Weight[i] = -1;
        }
    }
}

void bond2index(int L, int T) {
    if(Bond2index==NULL) {
        Bond2index = (int*)malloc(sizeof(int)*(T-1)*L*4);
    } else {
        free(Bond2index);
        Bond2index = (int*)malloc(sizeof(int)*(T-1)*L*4);
    }

    int ib=0;
    for(int t=1;t<T;t++) {
        for(int x=0;x<L;x++) {
            Bond2index[4*ib+0] = x+t*L;
            Bond2index[4*ib+1] = (x-1+L)%L + (t-1)*L;
            Bond2index[4*ib+2] = x + (t-1)*L;
            Bond2index[4*ib+3] = (x+1)%L + (t-1)*L;
            ib++;
        }
    }

    Nbond = ib;
    Nsite = L*T;
}

double WeightOfConf;
void clustering() {
    int s[4];
    int index[4];
    WeightOfConf=0;
    for(int ib=0;ib<Nbond;ib++) {
        index[0] = Bond2index[4*ib+0];
        index[1] = Bond2index[4*ib+1];
        index[2] = Bond2index[4*ib+2];
        index[3] = Bond2index[4*ib+3];

        s[0] = Sigma[index[0]];
        s[1] = Sigma[index[1]];
        s[2] = Sigma[index[2]];
        s[3] = Sigma[index[3]];

        if(gsl_rng_uniform_pos(rng)<Prob[bits2int(s)]) {
           merge(Ptree,Weight,index[0],index[1]); 
           merge(Ptree,Weight,index[0],index[2]); 
           merge(Ptree,Weight,index[0],index[3]); 
        }

        WeightOfConf += log(LocalWeight[bits2int(s)]);
    }
}

void flip() {
    for(int i=0;i<Nsite;i++) {
        Flist[i]=0;
        if(gsl_rng_uniform_pos(rng)<0.5) Flist[i]=1;
    }

    int r;
    for(int i=0;i<Nsite;i++) {
        r = root(Ptree,i);
        if(Weight[r]>0){
            Sigma[i] = (Sigma[i])^Flist[r];
        }
    }
}

int* Clusters;
void analysis_cluster() {
    if(Clusters==NULL) {
        Clusters = (int*)malloc(sizeof(int)*Nsite);
    }

    int r;
    int ncluster=0;
    int msize=0;
    for(int i=0;i<Nsite;i++) Clusters[i] = 0;
    for(int i=0;i<Nsite;i++) {
        r = root(Ptree,i);
        if(!Clusters[r]) {
            Clusters[r] = 1;
            ncluster++;
        } else {
            Clusters[r]++;

            if(Clusters[r]>msize) msize = Clusters[r];
        }
    }

    printf("--------------------------\n");
    printf("Nsite    = %d\n",Nsite);
    printf("ncluster = %d\n",ncluster);
    printf("msize    = %d\n",msize);
}

void show_state(int L, int T) {
    for(int t=0;t<T;t++) {
        for(int x=0;x<L;x++) {
            printf("%d ",Sigma[x+t*L]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    int L=atoi(argv[1]);
    int T=atoi(argv[2]);
    double beta=atof(argv[3]);
    int nblock  = atoi(argv[4]);
    int nsample = atoi(argv[5]);
    int nther   = atoi(argv[6]);
    unsigned long int seed = atoi(argv[7]);
    
    T+=1;

    Nsite = T*L;
    Nbond = (T-1)*L;

    Sigma = (int*)malloc(sizeof(int)*Nsite);
    Bond2index = (int*)malloc(sizeof(int)*Nbond*2);
    Ptree  = (int*)malloc(sizeof(int)*Nsite);
    Weight = (int*)malloc(sizeof(int)*Nsite);
    Flist  = (int*)malloc(sizeof(int)*Nsite);

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    bond2index(L,T);
    initial_state(L,T,2);

    int b[4];
    double w;
    for(int i=0;i<16;i++) {
        int2bits(i,b);
        w = configuration2weight(i,beta);
        LocalWeight[i] = w;
        printf("%d (%d | %d %d %d) w=%.6e\n",bits2int(b),b[0],b[1],b[2],b[3],w);
    }
    for(int i=0;i<16;i++) {
        Prob[i] = 1-LocalWeight[1]/LocalWeight[i];
        printf("%.6e \n",Prob[i]);
    }

    for(int i=0;i<nther;i++) {
        initial_tree();
        boundary_condition(L,T,0);
        clustering();
        flip();
        analysis_cluster();
    }

    printf("staring the measurement...\n");
    double* magz = (double*)malloc(sizeof(double)*T);

    char magz_filename[128];
    char weight_filename[128];
    sprintf(magz_filename,"magz_l_%d_beta_%.6f_t_%d.txt",L,beta,T);
    sprintf(weight_filename,"weight_l_%d_beta_%.6f_t_%d.txt",L,beta,T);
    FILE* magz_file = fopen(magz_filename,"w");
    FILE* weight_file = fopen(weight_filename,"w");

    for(int k=0;k<nblock;k++) {
        for(int i=0;i<T;i++) magz[i]=0;
        for(int i=0;i<nsample;i++) {
            for(int j=0;j<10;j++){
                initial_tree();
                boundary_condition(L,T,0);
                clustering();
                flip();
            }

            for(int t=0;t<T;t++) {
                for(int x=0;x<L;x++) {
                    magz[t]+=Sigma[x+L*t];
                }
            }
            if(i%10==0)
                fprintf(weight_file,"%16e ",WeightOfConf);
        }
        fprintf(weight_file,"\n");

        for(int t=0;t<T;t++) {
            magz[t] = magz[t]/nsample/L;
            magz[t] = magz[t]*2-1;

            fprintf(magz_file,"%.16e ",magz[t]);
        }
        fprintf(magz_file,"\n");
    }

    fclose(magz_file);
    fclose(weight_file);

    //for(int i=0;i<Nsite;i++) printf("%d %d\n",Ptree[i],Weight[i]);
    //show_state(L,T);

    return 0;
}
