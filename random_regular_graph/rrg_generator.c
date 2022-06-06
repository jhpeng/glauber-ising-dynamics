#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "union_find.h"

// random_regular_graph_generator 
//  This function generates random regular graph with total vetex number
//  n which has r edges for each vertex.
// 
//  parameter:
//   n (int) total vertex number
//   r (int) number of edge for each vertex
//   m (int) mode
//    - m=0 : no restriction
//    - m=1 : without self-looping
//    - m=2 : without self-looping a two or more edges between same 
//            vertices 
//
//   graph (int**) output of random regular graph
//   size (int*) output of the size of the graph
//   rng (gsl_rng) random number generator

int random_regular_graph_generator(int n, int r, int m, int* graph, gsl_rng* rng) {
    if(((r*n)%2)!=0) {
        printf("r*n should be an even number!\n");
        return 1;
    }

    int* nedge = (int*)malloc(sizeof(int)*n);

    int check_global=1;
    int ntrail=0;
    while(check_global) {
        for(int i=0;i<n;i++) nedge[i]=0;
        for(int i1=0;i1<n;i1++){
            for(int j1=nedge[i1];j1<r;){
                nedge[i1]++;

                int check_local=1;
                while(check_local){
                    int i2 = (int)(gsl_rng_uniform_pos(rng)*n);
                    //printf("i1 n1 i2 n2 : [%d %d %d %d]\n",i1,nedge[i1],i2,nedge[i2]);
                    if(nedge[i2]<r) {
                        int j2=nedge[i2];
                        nedge[i2]++;

                        graph[i1*r+j1] = i2;
                        graph[i2*r+j2] = i1;
                        check_local=0;
                    }
                }
                j1=nedge[i1];
            }
        }
        // checking for self-looping
        check_global=0;
        if(m==1 || m==2){
            for(int i=0;i<n;i++) {
                for(int j=0;j<r;j++) {
                    if(graph[i*r+j]==i) 
                        check_global=1;
                }
            }
        }
        // checking for two or more edges between same vertices
        if(m==2) {
            for(int i=0;i<n;i++) {
                for(int j1=0;j1<r;j1++){
                    for(int j2=j1+1;j2<r;j2++) {
                        if(graph[i*r+j1]==graph[i*r+j2])
                            check_global=1;
                    }
                }
            }
        }
        ntrail++;

        //printf("ntrail = %d",ntrail);
        //printf("\r");
        //fflush(stdout);
    }

    printf("ntrail = %d\n",ntrail);

    free(nedge);

    return 0;
}

int graph_analysis(int n, int r, int* graph) {
    int* p = (int*)malloc(sizeof(int)*n);
    int* w = (int*)malloc(sizeof(int)*n);

    for(int i=0;i<n;i++) {
        p[i]=i;
        w[i]=1;
    }

    for(int i=0;i<n;i++) {
        for(int j=0;j<r;j++) 
            merge(p,w,i,graph[i*r+j]);
    }

    int count=0;
    for(int i=0;i<n;i++) {
        if(p[i]==i) {
            w[count]=w[i];
            count++;
        }
    }

    printf("total number of isolated cluster : %d\n",count);
    printf("size of clusters [ ");
    for(int i=0;i<count;i++)
        printf("%d ",w[i]);

    printf("]\n");

    free(p);
    free(w);

    return count;
}
