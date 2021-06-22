
#include <math.h>
#include<stdio.h>
#include<time.h>
#include "mat_entries.h"
#include "states_creation.h"

#define IDX(kk,ll) kk*ll




int main(void)
{
    clock_t start = clock();
    int momentum=16;
    int n=113;  //dimension
    double coupling=1;
    int *pre_array;
    
    cudaMallocManaged(&pre_array,IDX(n,momentum)*sizeof(int)); //array initialization of partition of integer momentum
    
    for(int i=0;i<IDX(n,momentum);i++)
    {
	    pre_array[i]=0;
    }
    partition<<<1,1>>>(momentum,pre_array);
    cudaDeviceSynchronize();
    
   /*for(int i=0;i<IDX(n,momentum);i++)
    {
           printf("%i,",pre_array[i]);
    } */
    printf("\n");
    
    //main_array of making proper state 
    int *main_array;
    cudaMallocManaged(&main_array,IDX(n,28)*sizeof(int));   //every 28 places new state
    
    for(int i=0;i<IDX(n,28);i++)
    {
        main_array[i]=0; 
    }
    array_create<<<1,1>>>(pre_array ,main_array,n,momentum);
    cudaDeviceSynchronize();
   
    /*for(int i=0;i<n*28;i++)
    {
        
        if(main_array[i]!=0)
        {
            printf("%i,",main_array[i]);
        }
        
    }*/
  
    
    cudaFree(pre_array);//remove pre_array as we are done with it 
    
    double ii=(double) n/1024;
    printf("ii=%i\n",ii);
    ii=ceil(ii);
    printf("ii=%i\n",ii);
    double *mass_array;
    cudaMallocManaged(&mass_array,n*sizeof(double));
    mass<<<ii ,1024>>>(main_array,mass_array,n);
    cudaDeviceSynchronize();
    /*printf("mass elements");
    for(int i =0;i<n;i++)
    {
        printf("%f,",mass_array[i]);
    }*/
    double *diag_array;
    cudaMallocManaged(&diag_array,n*sizeof(double));
    diagonal<<<ii ,1024>>>(main_array,diag_array,n,coupling);
    cudaDeviceSynchronize();
    printf("diag elements");
    for(int i =0;i<n;i++)
    {
        printf("%f,",momentum*(mass_array[i]+diag_array[i]));
    }
    clock_t end = clock();
    printf("\ntime taken=%d",end - start);
    printf("\n");
    
     
    cudaFree(main_array);
    
    cudaFree(mass_array);
    
    
}
