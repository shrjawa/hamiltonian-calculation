#include "ndh.h"
#include <math.h>

double ndhamilt2(int* bc,int* comp ,double lambda,int num,int len,int *temparray)
{
    
    int i;
    int j;
    int g;
    int f;
    int t;
    int s;
    double comb1;
    double comb2;
    
    double NDH=0;
    for(i=0;i<len;i=i+2)
    {
        int k=bc[i];
        int kk=bc[i+1];
        if(k==0)
        {break;}
            for(j=i;j<len;j=j+2)
            {
                    
                    int l=bc[j];
                    
                    int ll=bc[j+1];
                    if(k==l)
                    {
                      comb1=1;
                    }
                    else{comb1=2;}

                    if(k!=0 && l!=0)
                    {

                        for(t=0;t<len;t=t+2)
                        {
                            int m=comp[t];
                            int mm=comp[t+1];
                            if(m==0){break;}
                            for(s=t;s<len;s=s+2)
                            {
                                    
                                    int n=comp[s];
                                    
                                    int nn=comp[s+1];
                                    if(n==0)
                                    {
                                        break;
                                    }
                                    if(m==n)
                                    {
                                      comb2=1;
                                    }
                                    else{comb2=2;}
                                    //int temparray[len+4]={0};
                                    if(k+l==m+n && m!=0 && n!=0)
                                    {
                                      
                                      //printf("chk03 value %i\n",chk03);
                                      int chk1=0;
                                      int chk2=0;
                                        for(f=0;f<len;f++)   //copying data into temparray which can be modified
                                        {
                                         temparray[f]=comp[f];
                                        }
                                        long double klmn=k*l*m*n;
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
                                        for(g=0;g<len;g=g+2)
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
                                            temparray[len]=k;
                                            temparray[len+1]++;
                                        }
                                        if(chk2==0)
                                        {
                                            temparray[len+2]=l;
                                            temparray[len+3]++;
                                        }
                                        for(g=0;g<len+4;g=g+2)
                                        {
                                            for(f=g+2;f<len+4;f=f+2)
                                            {
                                                if(temparray[g]<temparray[f])
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


                                      //  printf("\n");
                                        int yy=0;

                                        for(f=0;f<len;f++)
                                        {
                                            if(temparray[f]==bc[f])
                                            {
                                                yy++;
                                            }
                                        }
                                        
                                        if(yy==len && k!=0 && l!=0 && m!=0 && n!=0)
                                        {
                                            if(s==t)
                                            {
                                                nn--;
                                                NDH=num*lambda*(1/(16*M_PI))*(sqrt(mm))*(sqrt(nn))*(sqrt(kk))*(sqrt(ll))/(sqrt(klmn));
                                            }
                                            if(k==l)
                                            {
                                                ll--;
                                                NDH=num*lambda*(1/(16*M_PI))*(sqrt(mm))*(sqrt(nn))*(sqrt(kk))*(sqrt(ll))/(sqrt(klmn));
                                            }
                                            if(s!=t && k!=l)
                                            {
                                                NDH=num*lambda*(1/(16*M_PI))*(sqrt(mm))*(sqrt(nn))*(sqrt(kk))*(sqrt(ll))/(sqrt(klmn));
                                            }
                                            NDH=comb1*comb2*NDH;
                                            
                                        
                                        }
                                    }
                                if(NDH!=0)
                                {break;}
                            }
                            if(NDH!=0)
                            {break;}
                        }
                    }

            if(NDH!=0)
            {break;}
            }
        if(NDH!=0)
        {break;}
    }
 return NDH;   
}

double ndhamilt3(int *prin,int *creat,double lambda,int num, int len,int *temparray1)
{

    int i;
    int j;
    int g;
    int f;
    int t;
    int s;
    int comb1;
    
    double NDH=0;
    int chek=0;
    for (i=0;i<len;i=i+2)
    {
        int k=prin[i];
        int kk=prin[i+1];
        if(k==0)
        {
            break;
        }
        if(k!=0)
        {
            for(j=0;j<len;j=j+2)
            {
                int l=creat[j];
                int ll=creat[j+1];
                if(l==0)
                {
                    break;
                }
                for(g=j;g<len;g=g+2)
                {
                    int m=creat[g];
                    int mm=creat[g+1];
                    if(m==0)
                    {
                        break;
                    }
                    for(f=g;f<len;f=f+2)
                    {
                        
                        
                        
                        
                        int n=creat[f];
                        int nn=creat[f+1];
                        if(n==0)
                        {
                            break;
                        }
                        //printf("l=%i,ll=%i,m=%i,mm=%i,n=%i,nn=%i",l,ll,m,mm,n,nn);
                        if(l==m && m==n)
                        {
                          {comb1=1;}
                        }
                        if(l!=m && m==n){comb1=3;}
                        if(l==m && m!=n){comb1=3;}
                        if(l!=m && l==n){comb1=3;}
                        if(l!=m && m!=n){comb1=6;}
                       // int temparray1[len+2]={0};
                        int lmn=l+m+n;

                        if(l!=0 && m!=0 && n!=0)
                        {if(k==lmn)
                        {
                            for(s=0;s<len;s++)   //copying data into temparray which can be modified
                            {
                                temparray1[s]=creat[s];
                            }
                            chek++;
                            long double klmn1=k*l*m*n;
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
                            for(s=0;s<len;s=s+2)
                            {
                                if(k==temparray1[s])
                                {
                                    temparray1[s+1]++;
                                    chk1++;
                                }
                            }
                            if(chk1==0)
                            {
                                temparray1[len]=k;
                                temparray1[len+1]++;

                            }

                            for(s=0;s<len+2;s=s+2)
                            {
                                for(t=s+2;t<len+2;t=t+2)
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
                            for(t=0;t<len;t++)
                            {
                                if(temparray1[t]==prin[t])
                                {
                                    zz++;
                                }
                            }
                            long double repl=sqrt(ll)*sqrt(mm)*sqrt(mm-1)*sqrt(kk);
                            long double repl1=sqrt(klmn1);

                            if(zz==len)
                            {
                                if(g==f && f==j)
                                {
                                    NDH=num*lambda*(1/(24*M_PI))*(sqrt(ll))*(sqrt(ll-1))*(sqrt(ll-2))*(sqrt(kk))/(sqrt(klmn1));
                                   
                                }
                                if(j==g && j!=f)
                                {
                                    NDH=num*lambda*(1/(24*M_PI))*(sqrt(ll))*(sqrt(ll-1))*(sqrt(nn))*(sqrt(kk))/(sqrt(klmn1));
                                }
                                if(j==f && j!=g)
                                {
                                    NDH=num*lambda*(1/(24*M_PI))*(sqrt(ll))*(sqrt(ll-1))*(sqrt(mm))*(sqrt(kk))/(sqrt(klmn1));
                                }
                                if(g==f && g!=j)
                                {
                                    NDH=num*lambda*(1/(24*M_PI))*repl/repl1;
                                }
                                if(j!=g && g!=f && j!=f)
                                {
                                    NDH=num*lambda*(1/(24*M_PI))*(sqrt(ll))*(sqrt(mm))*(sqrt(nn))*(sqrt(kk))/(sqrt(klmn1));
                                }
                                NDH=NDH*comb1;
                                
                                
                                
                            }

                        }
                        }
                        if(NDH!=0)
                        {break;}
                    }
                if(NDH!=0)
                {break;}
                
                }
                if(NDH!=0)
                {break;}
            }
        }
        if(NDH!=0)
        {break;}

    }
return NDH;
}