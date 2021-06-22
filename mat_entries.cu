#include "mat_entries.h"

__global__
void mass(int *arr_m, double *mass_array,int n)
{
    int i = threadIdx.x;
    int k=blockIdx.x;
    int kk=k*1024;
    
    
    double dd=0;
    int oj,ii;
    if(i+kk<n)
    {
        ii=(i+kk)*28;
        for(oj=0;oj<28;oj=oj+2)
        {
            if(arr_m[ii+oj]!=0)
            {
                double mass_d;
                mass_d=(double) arr_m[ii+oj+1]/arr_m[ii+oj];
                dd=dd+mass_d;

            }
        }
        mass_array[i+kk]=dd;
    }
    //printf("mass=%f\n\n",dd);
} 

__global__
void diagonal(int *arr,double *diag_array,int n,double coupling)
{
    
    int j,ii;
    double a,s,c,v;
    
    int i = threadIdx.x;
    int k=blockIdx.x;
    int kk=k*1024;
    double diag_H=0;
    if(i+kk<n)    
    {
        ii=(i+kk)*28;
        for(j=0;j<28;j=j+2)
        {
            if(arr[ii+j]==0)
            {
                break;
            }
            if(arr[ii+j]!=0)
            {
                double dd1; 

                dd1=(double) (arr[ii+j+1]*(arr[ii+j+1]-1))/(arr[ii+j]*arr[ii+j]);
                diag_H=diag_H+ ((coupling*dd1)/(16*M_PI));
           
                for(k=j+2;k<28;k=k+2)
                {
                    if(arr[ii+k]!=0)
                    {
                        a=arr[ii+j];
                        s=arr[ii+k];
                        c=arr[ii+j+1];
                        v=arr[ii+k+1];
                        diag_H=diag_H+((4*coupling*c*v)/(a*s*16*M_PI));
                    }
                }

            }
        }
        diag_array[i+kk]=diag_H; 
    }
   // printf("diag_H=%f\n\n",diag_H);
   // return diag_H;
}
