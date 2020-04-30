#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

double num;
int len;
int dim=0;
int h=0;
void partitionAPBCodd(double n)
{

    double p[len];
    int  k = 0;
    p[k] = n;
    double temp;
    while (1)
    {
       /*for(int i =0;i<=k;i++)
	   {
		   printf("%f,",p[i]);
	   }*/
	   
         if(k%2==0)
			{
			for (int i=0; i<=k;i++)
			{
			int fch=2*p[i];
			if(fch%2==0)
			{h+=1;}
			}
			if (h==0){dim+=1;
			for(int i =0;i<=k;i++)
	   {
	   printf("%f,",p[i]);
	   if(i==k){printf("\n");}}
	   }
	   
			} 
		h=0;
        temp = 0;
        while (k >= 0 && p[k] == 0.5)
        {
            temp += 0.5;
            k--;
        }
		
        if (k < 0)  return;
		
        p[k]-=0.5;;
		//printf("%d"
        temp+=0.5;
		
        while (temp > p[k])
        {
            p[k+1] = p[k];
            temp = temp - p[k];
            k++;
        }
		//printf("temp=%f\n",temp);
		//printf("p[k]=%f\n",p[k]);
        p[k+1] = temp;
		//printf("p[k+1]=%f\n",p[k+1]);
        k++;
		//printf("k=%i\n",k);
    }

}


void partitionAPBCeven(double n)
{

    double p[len];
    int  k = 0;
    p[k] = n;
    double temp;
    while (1)
    {
       /*for(int i =0;i<=k;i++)
	   {
		   printf("%f,",p[i]);
	   }*/
	   
         if(k%2==1)
			{
			for (int i=0; i<=k;i++)
			{
			int fch=2*p[i];
			if(fch%2==0)
			{h+=1;}
			}
			if (h==0){dim+=1;
			for(int i =0;i<=k;i++)
	   {
	   printf("%f,",p[i]);
	   if(i==k){printf("\n");}}
	   }
	   
			} 
		h=0;
        temp = 0;
        while (k >= 0 && p[k] == 0.5)
        {
            temp += 0.5;
            k--;
        }
		
        if (k < 0)  return;
		
        p[k]-=0.5;;
		//printf("%d"
        temp+=0.5;
		
        while (temp > p[k])
        {
            p[k+1] = p[k];
            temp = temp - p[k];
            k++;
        }
		//printf("temp=%f\n",temp);
		//printf("p[k]=%f\n",p[k]);
        p[k+1] = temp;
		//printf("p[k+1]=%f\n",p[k+1]);
        k++;
		//printf("k=%i\n",k);
    }

}
int main()
{
    //printf("valu is %f\n",M_PI );
    printf("\nEnter a number to perform integer partition: ");
    scanf("%lf", &num);  //value of K.
    //printf("\nEnter renormalised mass value: ");
    //scanf("%lf",&mass);
    //num=16;
    //double mass=1.0;
	
	len=2*num;
	if(len%2==1)
    {partitionAPBCodd(num);
		}
	if(len%2==0)
    {partitionAPBCeven(num);
		}
		
	printf("%i\n",dim);
	}