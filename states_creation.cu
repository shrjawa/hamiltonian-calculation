#include "states_creation.h"

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


