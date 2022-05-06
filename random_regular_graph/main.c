#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "gsl/gsl_rng.h"

#include "rrg_generator.h"

int main() {
    int n=100;
    int r=6;
    unsigned long int seed=84893;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);

    int* graph;
    int size;

    if(random_regular_graph_generator(n,r,2,&graph,&size,rng))
        return 1;

    for(int i=0;i<n;i++) {
        printf("i=%d : [",i);
        for(int j=0;j<r;j++) {
            printf("%d ",graph[i*r+j]);
        }
        printf("]\n");
    }

    return 0;
}
