#include "stdio.h"
#include "stdlib.h"
#include "math.h"

static int stat_ndata;
static int stat_nbin;
static int stat_bsize;
static int stat_cap;
static int stat_len;
static int stat_index;
static double** stat_data;
static double** stat_data_bin;

void stat_setup(int ndata, int nbin) {
    stat_ndata  = ndata;
    stat_nbin   = nbin;
    stat_cap    = 10000000;
    stat_len    = 0;
    stat_index  = 0;
    stat_bsize  = stat_cap/nbin;

    stat_data = (double**)malloc(sizeof(double*)*ndata);
    stat_data_bin = (double**)malloc(sizeof(double*)*ndata);
    for(int i=0;i<ndata;i++) {
        stat_data[i] = (double*)malloc(sizeof(double)*stat_cap);
        stat_data_bin[i] = (double*)malloc(sizeof(double)*stat_nbin);
    }
}

void stat_append(double* data) {
    for(int i=0;i<stat_ndata;i++) {
        (stat_data[i])[stat_index] = data[i];
    }
    if(stat_len<stat_cap) {
        stat_len++;
    }

    stat_index++;
    if(stat_index==stat_cap) stat_index=0;
}

static double stat_ave(double* d, int n) {
    double ave=0;
    for(int i=0;i<n;i++) ave+=d[n];
    
    return ave/n;
}

static double stat_var(double* d, int n) {
    double ave = stat_ave(d,n);
    double var=0;
    for(int i=0;i<n;i++) var += (d[i]-ave)*(d[i]-ave);

    return var/n;
}

void stat_print_status() {
    double var[stat_ndata];
    double var_bin[stat_ndata];
    double tau[stat_ndata];
    int nbin=stat_len/stat_bsize;
    double* p;

    printf("---------------------------------\n");
    printf("status of Markov Chain\n");
    printf("stat_ndata = %d\n",stat_ndata);
    printf("stat_nbin / nbin = %d / %d\n",stat_nbin,nbin);
    printf("stat_bsize = %d\n",stat_bsize);
    printf("stat_cap = %d\n",stat_cap);
    printf("stat_len = %d\n",stat_len);
    if(nbin>1) {
    printf("             |     var    |  var_bin   |    tau\n");
    for(int i_data=0;i_data<stat_ndata;i_data++) {
        for(int i_bin=0;i_bin<nbin;i_bin++) {
            p = &((stat_data[i_data])[i_bin*stat_bsize]);
            (stat_data_bin[i_data])[i_bin] = stat_ave(p,stat_bsize);
        }

        var[i_data] = stat_var(stat_data[i_data],stat_len);
        var_bin[i_data] = stat_var(stat_data_bin[i_data],nbin);
        tau[i_data] = var_bin[i_data]*stat_bsize/var[i_data];

        printf("data id (%d) : %.6e %.6e %.6e\n",i_data,var[i_data],var_bin[i_data],tau[i_data]);
    }
    } else {
        printf("nbin should larger than one\n");
    }
    fflush(stdout);

}

void stat_free() {
    for(int i=0;i<stat_ndata;i++) {
        free(stat_data[i]);
        free(stat_data_bin[i]);
    }
    free(stat_data);
    free(stat_data_bin);
}
