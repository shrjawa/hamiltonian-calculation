#include "sort_arrays.h"
#include <stdbool.h>

void  bubblesort(int **main_array,int n,int len)
{
    int i,j,k;
    bool swpd;
    for(i=0;i<n;i++)
    {
        swpd=false;
        for(j=0;j<n-i-1;j++)
        {
            if(main_array[j][len]>main_array[j+1][len])
            {
                int temp;
                for(k=0;k<len+1;k++)
                {
                    temp=main_array[j][k];
                    main_array[j][k]=main_array[j+1][k];
                    main_array[j+1][k]=temp;
                }
                swpd=true;
            }
        }
        if(swpd==false)
        {
            break;
        }

    }
}