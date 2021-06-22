#ifndef states_creation_   /* Include guard */
#define states_creation_

__global__
void partition(int m,int *pre_array);

__global__
void array_create(int *pre_array,int *main_array,int n,int momentum);

#endif 