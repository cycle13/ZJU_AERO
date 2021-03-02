#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

float* get_refl(float *refl, int len, 
                float *Da, int len_D_bins_x,
                float *Db, int o1,
                float *rcs, int o2,
                float *N, int o3, 
                float step_D,
                float Dmin);

float* get_refl(float *refl, int len, 
                float *Da, int len_D_bins_x,
                float *Db, int o1,
                float *rcs, int o2,
                float *N, int o3, 
                float step_D,
                float Dmin){
   memset(refl, 0, len*sizeof(float));
   float sum=0;
   int i,j;
   int idx_D_a, idx_D_b;
   float D_a, D_b;
   
   for(i=0; i<len_D_bins_x; i++){
        sum=0;
        D_a=Da[i];
        D_b=Db[i];

        idx_D_a = (int)((D_a-Dmin)/step_D);
        idx_D_b = (int)((D_b-Dmin)/step_D);

        // printf("%d %d \n", idx_D_a, idx_D_b);

        for(j=idx_D_a; j<idx_D_b; j++){
            sum += N[j] * rcs[j];
        }
        refl[i] += sum * step_D;
    }
    return refl;
}
