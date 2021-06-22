
#include <math.h>
#include<stdio.h>
#include<time.h>
#include "mat_entries.h"
#define IDX(kk,ll) kk*ll



__global__
void partition(int m,int *pre_array)
{
    int j;
    int c=0;
    int *p=(int *)malloc(m);
    int ii=0;
    int  k = 0;
    p[k] = m;
    int temp;
    while (1)
    {
        
        if(k%2==0)
        {
            for(j=0;j<k+1;j++)
            {
                pre_array[j+c]=p[j];
            }
            ii++;
            c+=m;
        }
        temp = 0;
        while (k >= 0 && p[k] == 1)
        {
            temp += 1;
            k--;
        }

        if (k < 0)  break;

        p[k]--;
        temp++;

        while (temp > p[k])
        {
            p[k+1] = p[k];
            temp = temp - p[k];
            k++;
        }

        p[k+1] = temp;
        k++;
        
    }
    
}



__global__
void array_create(int *pre_array,int *main_array,int n,int momentum)
{
   
    int i;
    int j,k;
    int mai=0;
    
    for(i=0;i<n*momentum;i=i+momentum)
    {
        
        int kk=0;
        
        int previous_num=0;
        for(j=0;j<momentum;j++)
        {
            if(pre_array[i+j]!=0)
            {
                if(pre_array[i+j]==previous_num)
                {
                    main_array[mai+kk-1]+=1;
                }
                if(pre_array[i+j]!=previous_num)
                {
                    int appendornot=0;
                    int number=pre_array[i+j];
                    for(k=0;k<28;k=k+2)
                    {
                        if(number==main_array[mai+k])
                        {
                            appendornot+=1;
                        }
                        if(main_array[mai+k]==0)
                        {
                            break;
                        }
                    }
                    
                    if(appendornot==0)
                    {    
                        
                        main_array[mai+kk]=number;
                        main_array[mai+kk+1]=1;
                        kk+=2;
                    }
                    previous_num=pre_array[i+j];
                }
            }
            
        } 
       
        mai=mai+28;
        
    }
}



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
