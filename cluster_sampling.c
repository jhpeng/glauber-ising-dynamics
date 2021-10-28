#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "union_find.h"

int* Sigma;
int* Bond2index;
int Nsite,Nbond;
double P;

int* Ptree;
int* Weight;
int* Flist;

gsl_rng* rng;

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

void construct_bond2index(int L, int T) {
    int index1,index2;
    int ib=0;
    for(int t=1;t<T;t++) {
        for(int x=0;x<L;x++) {
            index1 = x+t*L;
            index2 = (x-1+L)%L + (t-1)*L;
            Bond2index[2*ib+0] = index1;
            Bond2index[2*ib+1] = index2;
            ib++;

            index2 = x + (t-1)*L;
            Bond2index[2*ib+0] = index1;
            Bond2index[2*ib+1] = index2;
            ib++;

            index2 = (x+1)%L + (t-1)*L;
            Bond2index[2*ib+0] = index1;
            Bond2index[2*ib+1] = index2;
            ib++;
        }
    }
/*
    printf("Nbond = %d\n",Nbond);
    for(int i=0;i<ib;i=i+3) {
        int i1 = Bond2index[2*i+0];
        int i2 = Bond2index[2*i+1];
        int i3 = Bond2index[2*i+2];
        int i4 = Bond2index[2*i+3];
        int i5 = Bond2index[2*i+4];
        int i6 = Bond2index[2*i+5];
        printf("(%d %d) (%d %d) (%d %d)\n",i1,i2,i3,i4,i5,i6);
    }
*/
}

void clustering() {
    int index1,index2;
    for(int ib=0;ib<Nbond;ib++) {
        index1 = Bond2index[2*ib];
        index2 = Bond2index[2*ib+1];

        if(Sigma[index1]==Sigma[index2]) {
            if(gsl_rng_uniform_pos(rng)<P) {
                merge(Ptree,Weight,index1,index2);
            }
        }
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
            if(Flist[r]) Sigma[i] = (Sigma[i])^1;
        }
    }
}

void show_state(int L, int T) {
    for(int t=0;t<T;t++) {
        for(int x=0;x<L;x++) {
            printf("%d ",Sigma[x+t*L]);
        }
        printf("\n");
    }
}

int main() {
    int L=32;
    int T=20;
    double beta=0.5;
    unsigned long int seed = 94732974;
    int nsample=10000;
    
    T+=1;
    P = 1-exp(-2*beta);

    Nsite = T*L;
    Nbond = (T-1)*L*3;

    Sigma = (int*)malloc(sizeof(int)*Nsite);
    Bond2index = (int*)malloc(sizeof(int)*Nbond*2);
    Ptree  = (int*)malloc(sizeof(int)*Nsite);
    Weight = (int*)malloc(sizeof(int)*Nsite);
    Flist  = (int*)malloc(sizeof(int)*Nsite);

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    construct_bond2index(L,T);
    initial_state(L,T,1);

    for(int i=0;i<nsample;i++) {
        initial_tree();
        boundary_condition(L,T,0);
        clustering();
        flip();
    }

    for(int i=0;i<Nsite;i++) printf("%d %d\n",Ptree[i],Weight[i]);
    show_state(L,T);
    printf("P=%.6e\n",P);

    return 0;
}
