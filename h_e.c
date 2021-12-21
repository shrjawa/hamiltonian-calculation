#include "h_e.h"
#include<math.h>

double mass(int *arr_m,int len)
{
    

    
    double dd=0;
    int oj;
    for(oj=0;oj<len;oj=oj+2)
    {
        if(arr_m[oj]!=0)
        {
            double mass_d;
            mass_d=(double) arr_m[oj+1]/arr_m[oj];
            dd=dd+mass_d;
            
        }
    }
   // printf("mass=%f\n\n",dd);
    return dd;
}    

double diagonal(int *arr,double coupling,int len)
{
    
    int j,k;
    double a,s,c,v;
    
    double diag_H=0;
    for(j=0;j<len;j=j+2)
    {
        if(arr[j]==0)
        {
            break;
        }
        if(arr[j]!=0)
        {
            double dd1; 
            
            dd1=(double) (arr[j+1]*(arr[j+1]-1))/(arr[j]*arr[j]);
            diag_H=diag_H+ ((coupling*dd1)/(16*M_PI));
        }
    }
    
    for(j=0;j<len;j=j+2)
    {
        if(arr[j]==0)
        {break;}
        if(arr[j]!=0)
        {
            for(k=j+2;k<len;k=k+2)
            {
                if(arr[k]!=0)
                {
                    a=arr[j];
                    s=arr[k];
                    c=arr[j+1];
                    v=arr[k+1];
                    diag_H=diag_H+((4*coupling*c*v)/(a*s*16*M_PI));
                }
            }
        
        }
    }
   // printf("diag_H=%f\n\n",diag_H);
    return diag_H;
}