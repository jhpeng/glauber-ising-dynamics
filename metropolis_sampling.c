#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int* Sigma;
double LocalWeight[16];

gsl_rng* rng;

int Nbond,Nsite;

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

int* Bond2index;
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

    printf("debug\n");
}

int* Index2bond;
void index2bond(int L, int T) {
    if(Bond2index==NULL) {
        Index2bond = (int*)malloc(sizeof(int)*T*L*4);
    } else {
        free(Index2bond);
        Index2bond = (int*)malloc(sizeof(int)*T*L*4);
    }
    
    for(int i=0;i<Nsite*4;i++) Index2bond[i] = -1;

    int i0,i1,i2,i3;
    for(int ib=0;ib<Nbond;ib++) {
        i0 = Bond2index[4*ib+0];
        i1 = Bond2index[4*ib+1];
        i2 = Bond2index[4*ib+2];
        i3 = Bond2index[4*ib+3];

        Index2bond[i0*4+0] = ib;
        Index2bond[i1*4+1] = ib;
        Index2bond[i2*4+2] = ib;
        Index2bond[i3*4+3] = ib;
    }
    printf("debug\n");
}

void localupdate(int L, int T) {
    int index = L+(int)(gsl_rng_uniform_pos(rng)*L*(T-1));
    int bond;
    int sigma[4];
    double wp,wn;
    double p=0;


    for(int i=0;i<4;i++) {
        bond = Index2bond[index*4+i];
        if(bond!=-1) {
            sigma[0] = Sigma[Bond2index[bond*4+0]];
            sigma[1] = Sigma[Bond2index[bond*4+1]];
            sigma[2] = Sigma[Bond2index[bond*4+2]];
            sigma[3] = Sigma[Bond2index[bond*4+3]];

            wp = LocalWeight[bits2int(sigma)];
            sigma[i] = (sigma[i])^1;
            wn = LocalWeight[bits2int(sigma)];

            p += (wn-wp);
        }
    }

    if(gsl_rng_uniform_pos(rng)<exp(p)) {
        Sigma[index] = (Sigma[index])^1;
    }
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
    int L = 32;
    int T = 41;
    double beta=0.5;
    int nsample = 100000;
    unsigned long int seed = 394934393454;

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    int i;
    int b[4];
    double w;
    for(i=0;i<16;i++) {
        int2bits(i,b);
        w = configuration2weight(i,beta);
        LocalWeight[i] = log(w);
        printf("%d (%d | %d %d %d) w=%.6e\n",bits2int(b),b[0],b[1],b[2],b[3],w);
    }

    bond2index(L,T);
    index2bond(L,T);

    Sigma = (int*)malloc(sizeof(int)*Nsite);
    initial_state(L,T,1);

    printf("===================== initial state =====================\n");
    show_state(L,T);

    for(i=0;i<nsample;i++) {
        for(int j=0;j<Nsite;j++)
            localupdate(L,T);
    }
    printf("======================= last state ======================\n");
    show_state(L,T);

    return 0;
}
