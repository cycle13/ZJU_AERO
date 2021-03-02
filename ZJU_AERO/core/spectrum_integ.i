 %module spectrum_integ
 %{
 /* Put header files here or function declarations like below */
extern float* get_refl(
    float *refl, int len, 
    float *Da, int len_D_bins_x, 
    float *Db, int o1,  
    float *rcs, int o2, 
    float *N, int o3, 
    float step_D, 
    float Dmin);
 %}
 
%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}


%apply (float* ARGOUT_ARRAY1, int DIM1) {(float *refl, int len)}
%apply (float* IN_ARRAY1, int DIM1) {(float *Da, int len_D_bins_x)}
%apply (float* IN_ARRAY1, int DIM1) {(float *Db, int o1)}
%apply (float* IN_ARRAY1, int DIM1) {(float *rcs,int o2)}
%apply (float* IN_ARRAY1, int DIM1) {(float *N,int o3)}


float* get_refl(
float *refl, int len, 
float *Da, int len_D_bins_x,
float *Db, int o1,
float *rcs, int o2,
float *N, int o3,
float step_D,
float Dmin);
