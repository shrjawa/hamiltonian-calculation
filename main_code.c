static char help[] = "Standard symmetric eigenproblem corresponding to the Laplacian operator in 1 dimension.\n\n""The command line options are:\n""  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <slepceps.h>
#include "ndh.h"
#include "h_e.h"
#include "sort_arrays.h"



double lambda;

void array_create(int *pre_array,int **main_array,int dim,int i,int ik,int size)
{
    len=size-1;
    int j,k;
    {
        int kk=0;
        
        int previous_num=0;
        for(j=0;j<ik;j++)
        {
            if(pre_array[j]!=0)
            {
                if(pre_array[j]==previous_num)
                {
                    main_array[i][kk-1]+=1;
                }
                if(pre_array[j]!=previous_num)
                {
                    int appendornot=0;
                    int number=pre_array[j];
                    for(k=0;k<len;k=k+2)
                    {
                        if(number==main_array[i][k])
                        {
                            appendornot+=1;
                        }
                        if(main_array[i][k]==0)
                        {
                            break;
                        }
                    }
                    
                    if(appendornot==0)
                    {    
                        main_array[i][kk]=number;
                        main_array[i][kk+1]=1;
                        kk+=2;
                    }
                    previous_num=pre_array[j];
                }
            }
        }
    }
    main_array[i][len]=ik;
}


void partition(int n,int **main_array,int dim,int size,int ooe)
{

    int j;
    int p[n];
    int ii=0;
    int  k = 0;
    p[k] = n;
    int temp;
    while (1)
    {
        
        if(k%2==ooe /*&& k<=7*/)
        { 
            { 
                array_create(p,main_array, dim,ii,k+1,size);
            }
            ii++;
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



main(int argc, char **argv)
{
  Mat            A; 
  EPS            eps;
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscInt       Istart,Iend,nev,maxit,its,nconv,n;
  PetscScalar    kr,ki;
  Vec            xr,xi;
  PetscErrorCode ierr;
  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);  
  if (ierr) return ierr;  //???
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  int num_procs, my_pid;
  int num_threads, thread_id;
  int i, j, k, p, q, col, index, check_nnz, local_nnz, global_nnz;  //note: global_nnz
  double start;
  //int momentum=28, n=1851;
  int momentum=10;
  n=20;
    int oddoreven=0;  //zero for odd sector;one for even
  //int momentum=48, n=73593;
  double lambda=36.2;
  its=20;
  tol=1e-12;
  maxit=250  ;
  nconv=3;
  nev=1;
  
  double tstart, tbefore, tafter, tend, time_spent;
    
  ierr = 0;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  tstart = MPI_Wtime();
  
  //printf("my_pid=%i\n",my_pid);
   ierr=MatCreate(PETSC_COMM_WORLD,&A);
   ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
   ierr = MatSetFromOptions(A);CHKERRQ(ierr);
   ierr = MatSetUp(A);CHKERRQ(ierr);  //check...

  if (my_pid==0)
    {
      printf("\n\n --- Next run --- \n\n");
      printf("lamda=%f\n",lambda);
      printf("momentum=%i\n",momentum);
      
      printf("number of MPI ranks=%i\n",num_procs);
      printf("\n start creating many-body basis\n");
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  tbefore = MPI_Wtime();
  /*create auxiliary array */
   //array initialization of partition of integer momentum

  int size_array=29;
  int *main_array[n]; //array initialization the way we woud like to use
  for (i=0;i<n;i++)
    {
      main_array[i]=(int*)calloc(size_array,sizeof(int));  //calloc allows to get zeros in places 28 is more than enough to reach K=90
    }
    partition(momentum,main_array,n,size_array,oddoreven);
    bubblesort(main_array,n,size_array-1);
   // array of states is created here
    for (i=0;i<n;i++)
    {
        for( j =0;j<29;j++)
        {
            if(main_array[i][j]!=0)
            {
                printf("%i,",main_array[i][j]);
            }
            
        }
        printf("\n");
    }
    
 int lkl=size_array+3;
    int *temp_array2;
    temp_array2=(int*) calloc(size_array+3,sizeof(int));   // new line
    int *temp_array3;
    temp_array3=(int*) calloc(size_array+1,sizeof(int));  //new line
  
  MPI_Barrier(MPI_COMM_WORLD);
  tafter = MPI_Wtime();  
  if (my_pid==0)
    {
      time_spent = tafter - tbefore;
      printf("time setup basis states  %f\n", time_spent);
      printf("\n start constructing matrix\n");
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  tbefore = MPI_Wtime();

 ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

  for(i=Istart;i<Iend;i++)
    {
        
        double entry=momentum * (mass(main_array[i],size_array-1)+diagonal(main_array[i],lambda,size_array-1));
        //printf("entry(%i)=%f\n",i,entry);
        ierr = MatSetValue(A,i,i,entry,INSERT_VALUES);CHKERRQ(ierr);  //Setting the diagonal entries.
    }



  local_nnz = 0;
  // diagonal matrix elements all nonzero
  double Diag[n];   
  for(i=0;i<n;i++)
    {
      if (i%num_procs != my_pid) continue;
      double entry=momentum * (mass(main_array[i],size_array-1)+diagonal(main_array[i],lambda,size_array-1));
      {
	Diag[i]=entry;
	local_nnz+=1;
      }
    }
  
  for(i=0;i<n;i++)  //comparing different row element in array,choose ith row
    {
      if (i%num_procs != my_pid) continue;
      p=main_array[i][28]; // p stores number of elements from array[i][]
      
      index=0;
      col;
      for(k=i+1;k<n;k++)  //chooses kTH row to compare with iTH row
	{
    
	  q=main_array[k][28]; //q stores number of elements from array[k][]
        if(q>p+2)
        {
            break;
        }
	  
	  if(p==q && p!=0)// && p==2)
	    {
	      double NDH2=ndhamilt2(main_array[i],main_array[k],lambda,momentum,size_array-1,temp_array2);//i is main array and k is comparee
	      if(NDH2!=0)
		{
              ierr = MatSetValue(A,i,k,NDH2,INSERT_VALUES);CHKERRQ(ierr);
              ierr = MatSetValue(A,k,i,NDH2,INSERT_VALUES);CHKERRQ(ierr);

		  index+=1;
		  local_nnz+=1;
		}
	    }
	  
	  if(q==p+2  && p!=0)
	    {
	      double NDH3=ndhamilt3(main_array[i],main_array[k],lambda,momentum,size_array-1,temp_array3);
	      if(NDH3!=0)
		{   
              ierr = MatSetValue(A,i,k,NDH3,INSERT_VALUES);CHKERRQ(ierr);
              ierr = MatSetValue(A,k,i,NDH3,INSERT_VALUES);CHKERRQ(ierr);


		  index+=1;
		  local_nnz+=1;
		}
	    }
	}


    }

  //time_spent = MPI_Wtime() - tbefore;
  //printf("local time constr.matrix %f\n", time_spent);
  //printf("local number of nonzeros %i\n", local_nnz);


  MPI_Reduce(&local_nnz, &global_nnz, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);  //To get nnz, use GLOBAL nnz.
  
  check_nnz = 2 * global_nnz + n;
  
  tafter = MPI_Wtime();  
  if (my_pid==0)
    {
      time_spent = tafter - tbefore;
      printf("number of MPI ranks      %i\n", num_procs);
      printf("time constructing matrix %f\n", time_spent);
      printf("total number of nonzeros %i\n\n", check_nnz);
    }

   
  
  for (i=0;i<n;i++)  //clear memory of allocated memory for the matrix creation.
    {
      free(main_array[i]);

    }  
  
  tend = MPI_Wtime();
  time_spent = tend - tstart;
  
 
  if (my_pid==0)
    {
      printf("momentum           = %i \n", momentum);
      printf("basis dimension    = %i \n", n);
      printf("number of nonzeros = %i \n", check_nnz);
      printf("\n Total matrix construction runtime %f on %i MPI ranks \n\n", time_spent, num_procs);
    }
  MPI_Barrier(MPI_COMM_WORLD);

  tstart=MPI_Wtime();
  if (my_pid==0)
    {
      printf("\n\n --- Next run --- \n\n");
      printf("lamda=%f\n",lambda);
      printf("momentum=%i\n",momentum);
      
      printf("number of MPI ranks=%i\n",num_procs);
      printf("\n start diagonalizing\n");
    }
     ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
     ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
     clock_t Start2=clock();//New
     //printf("matrix_assembly complete");
     ierr = MatCreateVecs(A,NULL,&xr);CHKERRQ(ierr);
     ierr = MatCreateVecs(A,NULL,&xi);CHKERRQ(ierr);
    
    ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
    
    ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
    
    ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
    
    ierr = EPSSolve(eps);CHKERRQ(ierr);
    //ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
    ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
    ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
    //ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
    
    ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);
    //nconv=3;
    ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);
    if (nconv>0) 
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"           k          ||Ax-kx||/||kx||\n""   ----------------- ------------------\n");CHKERRQ(ierr);
        for (i=0;i<nconv;i++) 
        {
            ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
            ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);
                #if defined(PETSC_USE_COMPLEX)
                re = PetscRealPart(kr);
                im = PetscImaginaryPart(kr);
                #else 
                re = kr;
                im = ki;
                #endif
                if (im!=0.0) 
                {
                    ierr = PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error); CHKERRQ(ierr);
                } 
                else 
                {
                    ierr = PetscPrintf(PETSC_COMM_WORLD,"   %0.16f       %12g\n",(double)re,(double)error); CHKERRQ(ierr);
                }
        }
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    }
    clock_t end2 = clock(); //New
    double time_spent2 = (double)(end2 - Start2) / CLOCKS_PER_SEC;  //This is the total time taken by the CPU.  New
    printf("time taken for diagonalization=%f\n",time_spent2);//New.
    //if (i%comm_size != rank)
    //{printf("time taken total=%f\n",time_spent2);}
    
    //time_t seconds1;
    //seconds1 = time(NULL);
    //printf("time taken time_t=%ld \n",(seconds1-seconds));
    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = VecDestroy(&xr);CHKERRQ(ierr);
    ierr = VecDestroy(&xi);CHKERRQ(ierr);
    
    clock_t end = clock();
    double time_spentTOT = (double)(end - start) / CLOCKS_PER_SEC; //New.
    printf("Total time taken=%f\n",time_spentTOT);//New.
    

tafter = MPI_Wtime();  
  if (my_pid==0)
    {
      time_spent = tafter - tbefore;
      printf("time spent diagonalizing  %f\n", time_spent);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  tend = MPI_Wtime();
  time_spent = tend - tstart;
  
 
  if (my_pid==0)
    {
      printf("momentum           = %i \n", momentum);
      printf("basis dimension    = %i \n", n);
      printf("number of nonzeros = %i \n", check_nnz);
      printf("\n Total runtime %f on %i MPI ranks \n\n", time_spent, num_procs);
    }
  //MPI_Finalize();
  ierr = SlepcFinalize();  //Is there a redundancy in this call?
    
  
  return ierr;

}
