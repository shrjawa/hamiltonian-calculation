#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int partition(int *arr,int size,int n)
{
    int len=0;
    int p[n];
    int  k = size-1;
    for(int i=0;i<size;i++)
    {
        p[i]=arr[i];
    }
    int temp;
    while (1)
    {
       
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
    
        len++;
    }
    return len+1;
}




int main(int argc, char **argv) {

   
    //define MPI 
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int num=atoi(argv[1]);
    int l_sum=0;
  
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i =0;i<num/2;i++)
    {
        if(i%size==rank)
        {
            if(i!=num/2-1)
            {
                int *arr=(int*) malloc(1*sizeof(int));
                arr[0]=2+i;
                int l=partition(arr,1,num);
                l_sum+=l;
                free(arr);   
            }
            if(i==num/2-1)
            {
                int *arr1=(int*) malloc(3*sizeof(int));
                arr1[0]=(num-2)/2;
                arr1[1]=(num-2)/2;
                arr1[2]=2;
                int l=partition(arr1,3,num);
                l_sum+=l;
                free(arr1);
            }
        } 
    }
     MPI_Barrier(MPI_COMM_WORLD);
    int global_sum;
    MPI_Allreduce(&l_sum, &global_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("Global sum: %i\n", global_sum+2);
    }
 
    MPI_Finalize();
    return 0;
}
