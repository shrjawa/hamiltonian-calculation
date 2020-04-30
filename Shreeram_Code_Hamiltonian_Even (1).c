//Compile: gcc -Wall -o /Users/mamoonsharaf/LF_Codes/Shreeram_Code_Hamiltonian_Even /Users/mamoonsharaf/LF_Codes/Shreeram_Code_Hamiltonian_Even.c  -lm -lgslcblas -lgsl 
//Run: ./Shreeram_Code_Hamiltonian_Even


//scp /Users/mamoonsharaf/Project/Project\ Files/Project1.c msharaf@hpc-class.its.iastate.edu:
//scp msharaf@hpc-class.its.iastate.edu:/home/msharaf/XmatDiag2525.mtx /Users/mamoonsharaf

//Nonzero matrix elements in line 655.
//Change line 713.
//Eigenvalues in line 741.
//Elements specified in line 790.
#include<stdio.h>  //issues in running...
#include<stdlib.h>
#include<math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
//#include "partition.c"
//int lib[15]={2,1,4,3,6,6,8,12,10,22,12,40,14,69,16,118,18,195,20,317,22,505,24,793,26,1224,28,1867,30,2811,32,4186,34,6168,36,9005,
//38,13026,40,18692,42,26613,44,37613,46,52815,48,73618,50,102162,52,140853,54,193144};
//int connection[22]={17, 11, 13, 7, 8, 9, 20, 12, 16, 15, 5, 4, 18, 3, 2, 19, 10, 1, 14, 0, 21, 6};
int o=0;
int w;
int num;
double mass;
int ro=0;
int co=0;
int array[1000000][50]={{0}};
int check=0;
int diagonalentrynumber=0;
int numhamilt=0;
long double dim=1224;
long double hamiltonian[1224][1224]={{0}};//make changes here
double coup=1.0 ;
int halfmomentum=0;
int halfmom=0;
int chk01=0;
int chk02=0;
int chk03=0;
double lambda=49.3453; 
//int uk=0;
int enr[13]={0}; //enr[K/2] make change here should be half the momentum value ie. num/2 as # of 2 particle states is always is num/2 (K/2);
int inr=0;
int unr=0;
int non_zero=0;
//int lib[15]={2,5,11,22,42,77,135,231,385,627,1002,1575,2436,3718,5604};


//2D array is made here from line 179

void writearray(int q[],int column)
{
    int i;
    for(i=0;i<column;i++)
    {
        array[ro][co]=q[i];
        //if(ro==11 || ro ==12 || ro==15){
            
         //   printf("%i,",array[ro][co]);}
        co++;
    }
    //printf("\n");
    //printf("\n");

    ro++;
    co=0;
}

// this part of programm counts nunmber of number in array and diagonal entries are created

void counting(int p[], int size,double mren)  //check: p[] is array of states? 
{
    double mass =mren;
    int count_occur(int p[], char map[], int num_elements, int start);
    void print_array(int p[], int num_elements);
    int numparticle=0;
    char *map=malloc(size *sizeof(int));
    halfmomentum++;
    double num_occ;
    int i;
    int q[20]={0};
    int b=0;
    double H;
    int r;
    int t;
    for (i = 0; i < size; i++)
    {

        if(i==0)
            {

                w=0;
            }
        w++;
        if(i==0)
            {
                H=0;
            }
        if (map[i] == 0)
        {

            b++;
            num_occ = count_occur(p, map, size, i);
            double u=p[i];
            double x=num_occ-1;
            H=H+( num_occ / u/mass )+(lambda*(1/(16*M_PI))*coup*num_occ*x/p[i]/p[i]);  //coup?
            q[(2*b)-2]=p[i];
            q[(2*b)-1]=num_occ;
            //printf("H1=%f\n",H);
        }
    }

    for(r=0;r<2*b;r=r+2)
            {
                for(t=r+2;t<2*b;)
                {
                double a=q[r];
                double s=q[t];
                double c=q[r+1];
                double v=q[t+1];
                H=H+(4*(lambda/(16*M_PI))*coup*c*v/a/s);
                t=t+2;
                //printf("H2=%f\n",H);
                }

            }
    diagonalentrynumber++;
    int g=2*b;
    //printf("g=%i,",g);
    writearray(q, g);
    H=num*H;
    hamiltonian[diagonalentrynumber-1][diagonalentrynumber-1]=H;
    //printf("hamiltonian[%i][%i]= %f\n",diagonalentrynumber-1,diagonalentrynumber-1,hamiltonian[diagonalentrynumber-1][diagonalentrynumber-1]);
    int halfnum=num/2;
    if(p[0]==halfnum && p[1]==halfnum)// p[2]==1 && p[3]==1)// && p[4]==2 && p[5]==2 &&p[6]==1&&p[7]==1 &&p[8]==1 && p[9]==1 && p[10]==1 && p[11]==1 && p[12]==1&& p[13]==1 && p[14]==1 &&p[15]==1&&p[16]==1&&p[17]==1)//&&p[18]==1&&p[19]==1&&p[20]==1)
    {
      halfmom=halfmomentum-1;
      //printf("halfmomentum array number is %i\n",halfmom );
    }
    if(size==2)
    {
      for(i=0;i<size;i++)
      {
        //printf("%i  ",p[i] );
      }
      //printf("\n");
      enr[inr]=halfmomentum-1;
      inr++;
    }

}

int count_occur(int p[], char map[], int num_elements, int start)  
/* checks array a for number of occurrances of value */
{
    int i, count = 0, value = p[start];

    for (i = start; i < num_elements; i++)
    {
        if (p[i] == value)
        {
            map[i] = 1;
            ++count; /* it was found */
        }
    }
    return (count);
}

// ends here

void partition(int n,double mren)
{

    int p[n];
    int  k = 0;
    p[k] = n;
    int temp;
    while (1)
    {
        double mass=mren;
        if(k%2==1)
        {counting(p,k+1,mass);}
        temp = 0;
        while (k >= 0 && p[k] == 1)
        {
            temp += 1;
            k--;
        }

        if (k < 0)  return;

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

//non diagonal contribution to hamiltonian due to  2nd term in hamiltonian

void ndhamilt2(int bc,int comp)
{
    //FILE*file;
    //file=fopen("hamiltonian.txt","a");
    int i;
    int j;
    int g;
    int f;
    int t;
    int s;
    double comb1;
    double comb2;
    double NDH=0;

    for(i=0;i<20;i=i+2)
    {
            for(j=i;j<20;j=j+2)
            {
                    int k=array[bc][i];
                    int l=array[bc][j];
                    int kk=array[bc][i+1];
                    int ll=array[bc][j+1];
                    if(k==l)
                    {
                      comb1=1;
                    }
                    else{comb1=2;}

                    if(k!=0 && l!=0)
                    {
                        for(t=0;t<20;t=t+2)
                        {
                            for(s=t;s<20;s=s+2)
                            {
                                    int m=array[comp][t];
                                    int n=array[comp][s];
                                    int mm=array[comp][t+1];
                                    int nn=array[comp][s+1];
                                    if(m==n)
                                    {
                                      comb2=1;
                                    }
                                    else{comb2=2;}
                                    int temparray[24]={0};
                                    if(k+l==m+n && m!=0 && n!=0)
                                    {
                                        int chk1=0;
                                        int chk2=0;
                                      chk03++;
                                      //printf("chk03 value %i\n",chk03);
                                        for(f=0;f<20;f++)   //copying data into temparray which can be modified
                                        {
                                         temparray[f]=array[comp][f];
                                        }
                                        double klmn=k*l*m*n;
                                        temparray[t+1]--;
                                        if(temparray[t+1]==0)
                                        {
                                            temparray[t]=0;
                                        }
                                        if(t==s && temparray[t+1]!=0)
                                        {
                                            temparray[t+1]--;
                                            if(temparray[t+1]==0)
                                            {
                                                temparray[t]=0;

                                            }
                                        }
                                        if(t!=s)
                                        {
                                            temparray[s+1]--;
                                        }
                                        if(temparray[s+1]==0)
                                        {
                                            temparray[s]=0;
                                        }
                                        for(g=0;g<20;g=g+2)
                                        {
                                            if(k==temparray[g])
                                            {
                                                temparray[g+1]++;
                                                chk1++;   
                                            }
                                            if(l==temparray[g])
                                            {
                                                temparray[g+1]++;
                                                chk2++;
                                            }
                                        }
                                        if(chk1==0)
                                        {
                                            temparray[20]=k;
                                            temparray[21]++;
                                        }
                                
                                        if(chk2==0)
                                        {
                                            temparray[22]=l;
                                            temparray[23]++;
                                        }
                                        
                                        for(g=0;g<24;g=g+2)
                                        {
                                            for(f=g+2;f<24;f=f+2)
                                            {
                                                if(temparray[g]<temparray[f])  //sorting the temparray by momentum
                                                {
                                                    int aa=temparray[g];
                                                    int bb=temparray[g+1];
                                                    temparray[g]=temparray[f];
                                                    temparray[g+1]=temparray[f+1];
                                                    temparray[f]=aa;
                                                    temparray[f+1]=bb;
                                                }
                                            }
                                        }
                                        
                                       /* if(bc==11 && comp==12){printf("temparray is \n");
                                            for (g=0;g<24;g++){
                                        printf("%i,",temparray[g]);
                                        }printf("i=%i,j=%i,s=%i,t=%i",i,j,s,t);}*/


                                       // printf("\n");
                                        int yy=0;
                                        

                                        for(f=0;f<20;f++)
                                        {
                                            if(temparray[f]==array[bc][f])
                                            {
                                                yy++;
                                            }
                                        }
                                        
                                       /* if(bc==11 && comp==12){
                                        printf("hamiltonian entry for  check %i\n",yy);}*/
                            
                                        if(yy==20 && k!=0 && l!=0 && m!=0 && n!=0)
                                        {
                                            if(s==t)
                                            {
                                                nn--;
                                                NDH=num*coup*lambda*(1/(16*M_PI))*(sqrt(mm))*(sqrt(nn))*(sqrt(kk))*(sqrt(ll))/(sqrt(klmn));  //kk: number of particles with momentum k.
                                            }
                                            if(k==l)
                                            {
                                                ll--;
                                                NDH=num*coup*lambda*(1/(16*M_PI))*(sqrt(mm))*(sqrt(nn))*(sqrt(kk))*(sqrt(ll))/(sqrt(klmn));
                                            }
                                            if(s!=t && k!=l)
                                            {
                                                NDH=num*coup*lambda*(1/(16*M_PI))*(sqrt(mm))*(sqrt(nn))*(sqrt(kk))*(sqrt(ll))/(sqrt(klmn));
                                            }
                                            //printf("hamiltonian(%i,%i) is %f\n",bc,comp,NDH);
                                            //printf("comb1 =%i,comb2=%i\n",comb1,comb2);
                                            hamiltonian[bc][comp]=comb1*comb2*NDH;
                                            hamiltonian[comp][bc]=comb2*comb1*NDH;
                                            numhamilt++;
                                        check--;
                                        }
                                    }
                            }
                        }
                    }


            }
    }
}

void ndhamilt3(int prin,int creat)
{

    int i;
    int j;
    int g;
    int f;
    int t;
    int s;
    double comb1;
    double comb2;
    double comb3;
    double NDH1=0;
    int chek=0;
    for (i=0;i<20;i=i+2)
    {
        int k=array[prin][i];
        int kk=array[prin][i+1];
        if(k!=0)
        {
            for(j=0;j<20;j=j+2)
            {
                for(g=j;g<20;g=g+2)
                {
                    for(f=g;f<20;f=f+2)
                    {
                        int l=array[creat][j];
                        int ll=array[creat][j+1];
                        int m=array[creat][g];
                        int mm=array[creat][g+1];
                        int n=array[creat][f];
                        int nn=array[creat][f+1];
                        //printf("l=%i,ll=%i,m=%i,mm=%i,n=%i,nn=%i",l,ll,m,mm,n,nn);
                        if(l==m && m==n)
                        {
                          {comb1=1;}
                        }
                        if(l!=m && m==n){comb1=3;}  //combinatoral factor chosen accordingly...
                        if(l==m && m!=n){comb1=3;}
                        if(l!=m && l==n){comb1=3;}
                        if(l!=m && m!=n){comb1=6;}
                        int temparray1[22]={0};
                        int lmn=l+m+n;

                        if(l!=0 && m!=0 && n!=0)
                        {if(k==lmn)
                        {
                            for(s=0;s<20;s++)   //copying data into temparray which can be modified
                            {
                                temparray1[s]=array[creat][s];
                            }
                            chek++;
                            double klmn1=k*l*m*n;
                            //printf("k=%i,l=%i,m=%i,n=%i,klmn1=%f\n",k,l,m,n,klmn1);
                            temparray1[j+1]--;
                            if(temparray1[j+1]==0)
                            {
                                temparray1[j]=0;
                            }
                            if(j==g && temparray1[j+1]!=0)
                            {
                                temparray1[j+1]--;
                                if(temparray1[j+1]==0)
                                {
                                temparray1[j]=0;
                                }
                                if(j==f && temparray1[j+1]!=0)
                                {
                                    temparray1[j+1]--;
                                    if(temparray1[j+1]==0)
                                    {
                                        temparray1[j]=0;
                                    }
                                }
                            }
                            if(j==f && temparray1[j+1]!=0)
                            {
                                if(j!=g)
                                {
                                  temparray1[j+1]--;
                                  if(temparray1[j+1]==0)
                                  {
                                    temparray1[j]=0;
                                  }
                                }
                            }
                            if(j!=g)
                            {
                                temparray1[g+1]--;
                                if(temparray1[g+1]==0)
                                {
                                    temparray1[g]=0;
                                }
                                if(g==f && temparray1[g+1]!=0)
                                {
                                    temparray1[g+1]--;
                                    if(temparray1[g+1]==0)
                                    {
                                        temparray1[g]=0;
                                    }
                                }
                            }
                            if(j!=f)
                            {
                                if(g!=f)
                                {
                                    temparray1[f+1]--;
                                    if(temparray1[f+1]==0)
                                    {
                                        temparray1[f]=0;
                                    }
                                }
                            }

                            int chk1=0;
                            for(s=0;s<20;s=s+2)
                            {
                                if(k==temparray1[s])
                                {
                                    temparray1[s+1]++;
                                    chk1++;
                                }
                            }
                            if(chk1==0)
                            {
                                temparray1[20]=k;
                                temparray1[21]++;

                            }

                            for(s=0;s<22;s=s+2)
                            {
                                for(t=s+2;t<22;t=t+2)
                                {
                                    if(temparray1[s]<temparray1[t])
                                    {
                                        int aa=temparray1[s];
                                        int bb=temparray1[s+1];
                                        temparray1[s]=temparray1[t];
                                        temparray1[s+1]=temparray1[t+1];
                                        temparray1[t]=aa;
                                        temparray1[t+1]=bb;
                                    }
                                }
                            }
                            int zz=0;
                            for(t=0;t<20;t++)
                            {
                                if(temparray1[t]==array[prin][t])
                                {
                                    zz++;
                                }
                            }
                            double repl=sqrt(ll)*sqrt(mm)*sqrt(mm-1)*sqrt(kk);
                            double repl1=sqrt(klmn1);

                            if(zz==20)
                            {
                                if(g==f && f==j)
                                {
                                    NDH1=num*coup*lambda*(1/(24*M_PI))*(sqrt(ll))*(sqrt(ll-1))*(sqrt(ll-2))*(sqrt(kk))/(sqrt(klmn1));
                                    //printf("ll=%i,kk=%i,NDH1=%f,coup=%f ,num=%i,klmn1=%f \n",ll,kk,NDH1,coup,num,klmn1);
                                    //printf("k=%i,l=%i,m=%i,n=%i\n",k,l,m,n);
                                }
                                if(j==g && j!=f)
                                {
                                    NDH1=num*coup*lambda*(1/(24*M_PI))*(sqrt(ll))*(sqrt(ll-1))*(sqrt(nn))*(sqrt(kk))/(sqrt(klmn1));
                                }
                                if(j==f && j!=g)
                                {
                                    NDH1=num*coup*lambda*(1/(24*M_PI))*(sqrt(ll))*(sqrt(ll-1))*(sqrt(mm))*(sqrt(kk))/(sqrt(klmn1));
                                }
                                if(g==f && g!=j)
                                {
                                    NDH1=num*coup*lambda*(1/(24*M_PI))*repl/repl1;
                                }
                                if(j!=g && g!=f && j!=f)
                                {
                                    NDH1=num*coup*lambda*(1/(24*M_PI))*(sqrt(ll))*(sqrt(mm))*(sqrt(nn))*(sqrt(kk))/(sqrt(klmn1));
                                }
                                // if(prin==0 && creat==3){
                                // printf("lambda=%f , repl=%f,repl1=%f\n",lambda,repl,repl1);
                                // printf("hamiltonian(%i,%i) is %f\n",prin,creat,NDH1);
                                // printf("comb1=%i\n",comb1);}
                                hamiltonian[prin][creat]=comb1*NDH1;
                                hamiltonian[creat][prin]=comb1*NDH1;
                                numhamilt++;
                                //  fprintf(file,"hamiltonian(%i,%i) is %f\n",prin,creat,NDH1);
                            }

                        }}
                    }
                }
            }
        }

    }

}

int main()
{
    //printf("valu is %f\n",M_PI );
    printf("Enter a number to perform integer partition: ");
    scanf("%d", &num);  //value of K.
    double mass=1;
   // printf("\nEnter renormalised mass value: ");  //value of m.
   // scanf("%lf",&mass);  //long float...

    partition(num,mass);  //produces the states: at array[][]?
    //uk=num/2;
    //int *enr=malloc(sizeof(int)*uk);
    int i;
    int j;
    int k;
    //int *coef=malloc(num/2 *sizeof(int));


    for(i=0;i<dim;i++)  //comparing different row element in array,choose ith row
    {
        int p=0; // p stores number of elements from array[i][]    //is that the array of states?
        if(array[i][0]!=0)   
        {
            for(j=1;j<20;j=j+2)
            {
                p=p+array[i][j]; //???
            }
            //if(p==2)
            //{
              //coef[i]=i;
              //i++;
            //}
        }

        int col;
        
        for(k=i+1;k<dim;k++)  //chooses kTH row to compare with iTH row.  dim is defined at the very start to be 69.
        {
            int q=0; //q stores number of elements from array[k][]
            if(array[k][0]!=0)
            {
                for(col=1;col<20;col=col+2)
                {
                    q=q+array[k][col];  //???
                }
            }
            if(p==q && p!=0)// && p==2)
            {
                chk01++;   //chk01??
                ndhamilt2(i,k);//i is main array and k is compared.  How it functions??

            }

            if(q==p+2  && p!=0)
            {
                chk02++;
                ndhamilt3(i,k);
            }
        }
    }

   // printf("\n" );
    for(i=0;i<dim;i++)
    {
      for(j=i;j<dim;j++)
      {
        if (hamiltonian[i][j]!=0) {
          /* code */
            non_zero=non_zero+1;
      }
      //printf("\n");
    }}
    int rowCount=0;


    /*printf("Nonzero Matrix elements are: \n");
    for(i=0;i<dim;i++)
    {
        for(j=0;j<=i;j++)
        {
        if(hamiltonian[i][j]!=0)
        {
        printf("%f\n",hamiltonian[i][j]/8);
        }
        rowCount+=1;
        }
    
    }*/

    //printf("The number of rows is: %i, \n",rowCount);
    //double non_z=non_zero/3/3;
    printf("The number of nonzero elements in the upper triangular matrix are:\n");
    printf("%i",non_zero);
    //printf("%f,%f\n",non_zero,non_z);
   // printf("\n");
 gsl_matrix *m = gsl_matrix_alloc (dim,dim);//make changes here
 for (i=0;i<dim;i++)
  {
     for(j=0;j<dim;j++)
      {
           gsl_matrix_set (m, i, j, hamiltonian[i][j]);
      }
  }

//Write in the eigenvectors into input files.


/*FILE *fpVec; 

fpVec=fopen("mSqDiagK4K4.mtx","w"); 

for(i=0; i<dim; i++)
{ 
gsl_vector_view m_i= gsl_matrix_column (m, i);
for(j=0; j<dim; j++)
{
fprintf(fpVec, " %g ", gsl_vector_get(&m_i.vector,j)); /*use gsl data type*/
//}
 //fprintf(fpVec, "\n");
//}
//fclose(fpVec); 

//Write Diagonal Entries 

//FILE *fpVecOne; 

//fpVecOne=fopen("mSqDiagK4K4.mtx","w"); 

/*for(i=0; i<dim; i++)
{ 


fprintf(fpVecOne, "%f,",hamiltonian[i][i]/4);

}
fclose(fpVecOne); */

     gsl_vector *eval = gsl_vector_alloc (dim);//make changes here
     gsl_matrix *evec = gsl_matrix_alloc (dim,dim);

     gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (dim);

     gsl_eigen_symmv (m, eval, evec, w);  //I see that you are using the gsl symmetric diagonalization techniques here.

     gsl_eigen_symmv_free (w);
     gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
     double getmax=gsl_vector_get(eval,0);
     double maxelem=0;
     for (i=0;i<dim;i++)
        {
          gsl_vector *v=gsl_vector_alloc(dim);
          gsl_matrix_get_col (v, evec, i);
          float maxele=gsl_vector_get(v,halfmom);
          float maxele1=fabsf(maxele);
          //printf("%f\n",maxele1);
          if(maxele1>maxelem)
          {
            maxelem=maxele1;
          }
        }
    printf("The eigenvalues are: \n");
      for (i=0;i<10;i++)
      {
        double eigen_v=gsl_vector_get(eval,i);
        printf("%.16g\n",eigen_v); //***Prints eigenvalues.
        gsl_vector *ui=gsl_vector_alloc(dim);
        gsl_matrix_get_col(ui,evec,i);
        //printf("eigenvector is");
       // for (int j=0;j<dim;j++)
       // {
          //printf(" %g  ",gsl_vector_get(ui,j));
       // }
        //printf("\n\n\n\n");
      }
    printf("\n");
    printf("The lowest eigenvalue for lambda= %f is: %.16g\n ",lambda,gsl_vector_get(eval,0));

    // {
       //int i;
       //int j;
       double c1=0;
       double c3=0;
      /* for (i = 0; i < dim; i++)
       {
         double eval_i = gsl_vector_get (eval, i);
         gsl_vector *v=gsl_vector_alloc(dim);
         gsl_matrix_get_col (v, evec, i);
         double element = gsl_vector_get(v,halfmom);
         double element1=fabsf(element);
        /* if(maxelem==element1)
         {
         {
           //printf ("desired eigenvalue without normalization =%g\n", eval_i);
           for (j=0;j<num/2;j++)
           {

             //printf("evect(%i) %g\n",j,gsl_vector_get(v,enr[j]) );
             c1=c1+((gsl_vector_get(v,enr[j]))*(gsl_vector_get(v,enr[j])));// adding coeficient square of 2 particle states only
             //printf("value of c1=%.16f\n",c1);
           }
           for(j=0;j<dim;j++)
           {
             c3=c3+((gsl_vector_get(v,j))*(gsl_vector_get(v,j)));  // this should always yield 1 as it is checking orthonormalization condition
           }
        }
        // printf("eigenvector elements are\n");
        //  for(j=0;j<dim;j++)
        //   {
        //
        //     printf("%g\n",gsl_vector_get(v,j));
        //   }
        // printf("value of c0=%.16f\n",maxelem*maxelem);
        // printf("value of c1=%.16f\n",c1);
        // printf("c2 is %.16f",c3);  //prints c1 & c2.
        // }
     // }
    //  }
      //printf("\n" );
      /*printf("The nonzero Hamiltonian elements are: \n");
        for (i=0;i<dim;i++)
        {
            int con= connection[i];
          for (j=i;j<dim;j++)
          {
              int nnect=connection[j];
              double haml=hamiltonian[con][nnect]/4;
              if(haml!=0)
              {
                  printf("%i,%i , ",i+1,j+1);
            printf("%f , \n ",haml);
              }
          }
         printf("\n");
        }*/

      gsl_vector_free (eval);
      gsl_matrix_free (evec);
      // for (i=0;i<1000000;i++)
      //   {
      //     for(j=0;j<50;j++)
      //     {array[i][j]=0;}
      //   }
       // for(i=0;i<dim;i++)
       //
       // {
       //   for(j=0;j<i;j++)
       //   {
       //    if(hamiltonian[i][j]==0)
       //    {
       //      printf("N  ");
       //    }
       //    else
       //    {
       //      printf("Y  ");
       //    }
       //
       //   }
       //   printf("\n");
       // }
      // diagonalentrynumber=0;
    //  mass=getmax;
      //printf("mass renormalization value %lf\n",mass);
      // partition(num,mass);
      // for(i=0;i<dim;i++)  //comparing different row element in array,choose ith row
      // {
      //     int p=0; // p stores number of elements from array[i][]
      //     if(array[i][0]!=0)
      //     {
      //         for(j=1;j<20;j=j+2)
      //         {
      //             p=p+array[i][j];
      //         }
      //
      //     }
      //
      //     int col;
      //     for(k=i+1;k<dim;k++)  //chooses kTH row to compare with iTH row
      //     {
      //         int q=0; //q stores number of eelements from array[k][]
      //         if(array[k][0]!=0)
      //         {
      //             for(col=1;col<20;col=col+2)
      //             {
      //                 q=q+array[k][col];
      //             }
      //         }
      //         if(p==q && p!=0)// && p==2)
      //         {
      //             chk01++;
      //             ndhamilt2(i,k);//i is main array and k is comparee
      //             printf("value of p %i\n",p);
      //         }
      //
      //         if(q==p+2  && p!=0)
      //         {
      //             chk02++;
      //             ndhamilt3(i,k);
      //         }
      //     }
      // }
      // for (i=0;i<1000000;i++)
      //    {
      //        for(j=0;j<20;j++)
      //        {
      //            array[i][j]=0;
      //        }
      //    }
      // printf("checking hamiltonian");
      // for(i=0;i<dim;i++)
      // {
      //   for(j=0;j<dim;j++)
      //   {
      //     printf("%f,",hamiltonian[i][j]);
      //   }
      //   printf("\n");
      // }
   // gsl_matrix *m1 = gsl_matrix_alloc (dim,dim);//make changes here
   // for (i=0;i<dim;i++)
   //  {
   //     for(j=0;j<dim;j++)
   //      {
      //        gsl_matrix_set (m1, i, j, hamiltonian[i][j]);
   //      }
   //  }
      //
      //
      //  gsl_vector *eval1 = gsl_vector_alloc (dim);//make changes here
      //  gsl_matrix *evec1 = gsl_matrix_alloc (dim,dim);
      //
      //  gsl_eigen_symmv_workspace * w1 = gsl_eigen_symmv_alloc (dim);
      //
      //  gsl_eigen_symmv (m1, eval1, evec1, w1);
      //
      //  gsl_eigen_symmv_free (w1);
      //  gsl_eigen_symmv_sort (eval1, evec1, GSL_EIGEN_SORT_ABS_ASC);
      //  double getmax1=gsl_vector_get(eval1,0);
      //  printf("%f\n",getmax1);
      //  double maxelem1=0;
      //  for (i=0;i<dim;i++)
      //     {
      //       gsl_vector *v1=gsl_vector_alloc(dim);
      //       gsl_matrix_get_col (v1, evec1, i);
      //       float maxele2=gsl_vector_get(v1,1);
      //       float maxele3=fabsf(maxele2);
      //       //printf("maxele is %f\n",maxele3);
      //       if(maxele3>maxelem1)
      //       {
      //         maxelem1=maxele3;
      //       }
      //     }
      //
      //
      //  {
      //    int i;
      //    for (i = 0; i < dim; i++)
      //    {
      //      double eval1_i = gsl_vector_get (eval1, i);
      //      gsl_vector *v1=gsl_vector_alloc(dim);
      //      gsl_matrix_get_col (v1, evec1, i);
      //      double element2 = gsl_vector_get(v1,1);
      //      double element3=fabsf(element2);
      //      if(maxelem1==element3)
      //      {
      //      {
      //        printf ("desired eigenvalue =%g\n", eval1_i);
      //      }
      //      }
      //      //gsl_vector_fprintf (stdout, v, "%g");
      //    }
      //   }


      //gsl_vector_free (eval1);
      //gsl_matrix_free (evec1);

     return 0;
}


