#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "union_find.h"
#include "rrg_generator.h"

int bits2int(int nleg, int *b) {
    int i=0;
    for(int j=0;j<nleg;j++) {
        i+=(b[j])<<j;
    }

    return i;
}

void int2bits(int i, int nleg, int* b) {
    for(int j=0;j<nleg;j++) {
        b[j] = i&1;
        i = i>>1;
    }
}

double configuration2weight(int i, int nleg, double beta) {
    int b[nleg],s;

    int2bits(i,nleg,b);

    s=0;
    for(int j=1;j<nleg;j++) s+=(2*b[j]-1);
    double w = exp((b[0]*2-1)*s*beta);

    return w/(w+1.0/w);
}

void probability_adding_graph(int nleg, double beta, double* probs) {
    int b[nleg];
    int k=1<<nleg;
    double w;
    double w0 = configuration2weight(1,nleg,beta);

    for(int i=0;i<k;i++) {
        w = configuration2weight(i,nleg,beta);
        probs[i] = 1.0 - w0/w;

        int2bits(i,nleg,b);
        printf("%d (%d|",i,b[0]);
        for(int j=1;j<nleg;j++) printf(" %d",b[j]);
        printf(") w=%.10e\n",w);
    }
}

void initialization_tree(int nsite, int* p, int* w) {
    for(int i=0;i<nsite;i++) {
        p[i]=i;
        w[i]=1;
    }
}

void clustering(int nbond, int nleg, int* state, int* p, int* w, int* bond2index, double* probs, gsl_rng* rng) {
    int s[nleg];
    int index[nleg];

    for(int ib=0;ib<nbond;ib++) {
        for(int i=0;i<nleg;i++) index[i]=bond2index[nleg*ib+i];
        for(int i=0;i<nleg;i++) s[i]=state[index[i]];

        if(gsl_rng_uniform_pos(rng)<probs[bits2int(nleg,s)]) {
            for(int i=1;i<nleg;i++) 
                merge(p,w,index[0],index[i]);
        }
    }
}

void flipping(int nsite, int* state, int* p, int* w, int* flist, gsl_rng* rng) {
    for(int i=0;i<nsite;i++) {
        flist[i]=0;
        if(gsl_rng_uniform_int(rng,2)) flist[i]=1;
    }

    int r;
    for(int i=0;i<nsite;i++) {
        r=root(p,i);
        if(w[r]>0) {
            state[i] = (state[i])^flist[r];
        }
    }
}

void create_bond2index(int n, int r, int t_max, double beta, int* nsite, int* nbond, int* bond2index, double* probs, gsl_rng* rng) {
    int nleg=r+2;
    int* graph=(int*)malloc(sizeof(int)*n*r);

    printf("dtsw_method.c 98\n");
    random_regular_graph_generator(n,r,2,graph,rng);
    while(graph_analysis(n,r,graph)!=1) {
        random_regular_graph_generator(n,r,2,graph,rng);
    }

    int ib=0;
    for(int t=1;t<(t_max+1);t++) {
        for(int i=0;i<n;i++) {
            bond2index[nleg*ib+0] = i+t*n;
            bond2index[nleg*ib+1] = i+(t-1)*n;
            for(int j=0;j<r;j++)
                bond2index[nleg*ib+j+2] = graph[i*r+j]+(t-1)*n;

            ib++;
        }
    }

    *nbond = ib;
    *nsite = n*(t_max+1);

    probability_adding_graph(nleg,beta,probs);

    free(graph);
}

void initial_state(int n, int t_max, int type, int* state, gsl_rng* rng) {
    int nsite = n*(t_max+1);
    for(int i=0;i<nsite;i++) {
        state[i]=0;
        if(gsl_rng_uniform_int(rng,2)) state[i]=1;
    }
    for(int i=(n*t_max);i<nsite;i++)
        state[i]=1;

    if(type==1) {
        for(int i=0;i<n;i++)
            state[i]=0;
    }
}

void boundary_condition(int n, int t_max, int type, int nfix, double p, int* state, int* w, gsl_rng* rng){
    if(type==0) {
        for(int i=0;i<n;i++) {
            if(state[i]==1 && gsl_rng_uniform_pos(rng)<(2.0-1.0/p)) w[i]=-1;
        }
        for(int i=0;i<nfix;i++) {
            if(state[t_max*n+i]==1) 
                w[t_max*n+i]=-1;
        }
    } else if(type==1) {
        for(int i=0;i<n;i++) {
            if(state[i]==0) w[i]=-1;
        }
        for(int i=0;i<n;i++) {
            if(state[t_max*n+i]==1) 
                w[t_max*n+i]=-1;
        }
    }
}

void dtsw_cluster_analysis(int n, int t_max, int* p, int* w) {
    int ncluster=0;
    printf("cluster analysis\n");
    printf("cluster size : [");
    for(int i=0;i<(t_max+1)*n;i++) {
        if(i==p[i]) {
            ncluster++;
            printf("%d ",w[i]);
        }
    }
    printf("]\n");
    printf("ncluster = %d\n",ncluster);
    printf("------------------------------\n");
}

// only be used in the following function
static int dtsw_n, dtsw_r, dtsw_t_max, dtsw_type, dtsw_nfix,dtsw_nsite, dtsw_nbond, dtsw_nleg;
static int* dtsw_bond2index=NULL;
static int* dtsw_ptree=NULL;
static int* dtsw_weight=NULL;
static int* dtsw_state=NULL;
static int* dtsw_flist=NULL;
static double* dtsw_probs=NULL;
static double dtsw_beta,dtsw_p;


// measurement
static int dtsw_count=0;
static double* dtsw_measure_mz=NULL;
static double* dtsw_measure_mz_abs=NULL;
static double* dtsw_measure_mz_1_one=NULL;
static double* dtsw_measure_mz_1_two=NULL;

void dtsw_setup(int n, int r, int t_max, int type, int nfix, double beta, double p, gsl_rng* rng) {
    // parameters
    dtsw_n = n;
    dtsw_r = r;
    dtsw_nleg=r+2;
    dtsw_t_max = t_max;
    dtsw_type = type;
    dtsw_nfix = nfix;
    dtsw_beta = beta;
    dtsw_p = p;

    int nleg=r+2;
    dtsw_bond2index = (int*)malloc(sizeof(int)*nleg*n*t_max);
    dtsw_probs = (double*)malloc(sizeof(double)*nleg);
    create_bond2index(dtsw_n,dtsw_r,dtsw_t_max,dtsw_beta,&dtsw_nsite,&dtsw_nbond,dtsw_bond2index,dtsw_probs,rng);

    // placeholder
    dtsw_ptree = (int*)malloc(sizeof(int)*dtsw_nsite);
    dtsw_weight = (int*)malloc(sizeof(int)*dtsw_nsite);
    dtsw_state = (int*)malloc(sizeof(int)*dtsw_nsite);
    dtsw_flist = (int*)malloc(sizeof(int)*dtsw_nsite);

    initial_state(dtsw_n,dtsw_t_max,type,dtsw_state,rng);

    // measurement
    dtsw_measure_mz = (double*)malloc(sizeof(double)*(t_max+1));
    dtsw_measure_mz_abs = (double*)malloc(sizeof(double)*(t_max+1));
    dtsw_measure_mz_1_one = (double*)malloc(sizeof(double)*(t_max+1));
    dtsw_measure_mz_1_two = (double*)malloc(sizeof(double)*(t_max+1)*(t_max+1));
    for(int t=0;t<(dtsw_t_max+1);t++) {
        dtsw_measure_mz[t]=0;
        dtsw_measure_mz_abs[t]=0;
        dtsw_measure_mz_1_one[t]=0;
        for(int s=0;s<(dtsw_t_max+1);s++) {
            dtsw_measure_mz_1_two[t*(dtsw_t_max+1)+s]=0;
        }
    }
}

void dtsw_update(gsl_rng* rng) {
    initialization_tree(dtsw_nsite,dtsw_ptree,dtsw_weight);
    boundary_condition(dtsw_n,dtsw_t_max,dtsw_type,dtsw_nfix,dtsw_p,dtsw_state,dtsw_weight,rng);
    clustering(dtsw_nbond,dtsw_nleg,dtsw_state,dtsw_ptree,dtsw_weight,dtsw_bond2index,dtsw_probs,rng);
    flipping(dtsw_nsite,dtsw_state,dtsw_ptree,dtsw_weight,dtsw_flist,rng);
}

void dtsw_measurement_sampling() {
    double mz;
    for(int t=0;t<(dtsw_t_max+1);t++) {
        mz=0;
        for(int i=0;i<dtsw_n;i++) {
            mz+=dtsw_state[t*dtsw_n+i];
        }
        mz = (mz*2-dtsw_n)/dtsw_n;

        dtsw_measure_mz[t] += mz;
        dtsw_measure_mz_abs[t] += fabs(mz);
        dtsw_measure_mz_1_one[t] += dtsw_state[t*dtsw_n+0]*2-1;
        for(int s=0;s<(dtsw_t_max+1);s++) {
            double mz_t = dtsw_state[t*dtsw_n+0]*2-1;
            double mz_s = dtsw_state[s*dtsw_n+0]*2-1;
            dtsw_measure_mz_1_two[t*(dtsw_t_max+1)+s] += mz_t*mz_s;
        }
    }

    dtsw_count++;

    if(dtsw_count%1000==0) {
        printf("dtsw_count : %d \n",dtsw_count);
        dtsw_cluster_analysis(dtsw_n,dtsw_t_max,dtsw_ptree,dtsw_weight);
    }
}

void dtsw_measurement_save() {
    FILE* mz_file = fopen("dtsw_mz.txt","a");
    FILE* mz_abs_file = fopen("dtsw_mz_abs.txt","a");
    FILE* mz_1_one_file = fopen("dtsw_mz_1_one.txt","a");
    FILE* mz_1_two_file = fopen("dtsw_mz_1_two.txt","a");

    for(int t=0;t<(dtsw_t_max+1);t++) {
        fprintf(mz_file,"%.16e ", dtsw_measure_mz[t]/dtsw_count);
        fprintf(mz_abs_file,"%.16e ", dtsw_measure_mz_abs[t]/dtsw_count);
        fprintf(mz_1_one_file,"%.16e ", dtsw_measure_mz_1_one[t]/dtsw_count);

        dtsw_measure_mz[t]=0;
        dtsw_measure_mz_abs[t]=0;
        dtsw_measure_mz_1_one[t]=0;

        for(int s=0;s<(dtsw_t_max+1);s++) {
            fprintf(mz_1_two_file,"%.16e ", dtsw_measure_mz_1_two[t*(dtsw_t_max+1)+s]/dtsw_count);

            dtsw_measure_mz_1_two[t*(dtsw_t_max+1)+s]=0;
        }
    }
    fprintf(mz_file,"\n");
    fprintf(mz_abs_file,"\n");
    fprintf(mz_1_one_file,"\n");
    fprintf(mz_1_two_file,"\n");

    fclose(mz_file);
    fclose(mz_abs_file);
    fclose(mz_1_one_file);
    fclose(mz_1_two_file);

    dtsw_count=0;
}
