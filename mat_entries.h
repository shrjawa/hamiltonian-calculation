#ifndef mat_entries_   /* Include guard */
#define mat_entries_

__global__
void mass(int *arr_m, double *mass_array,int n);  /* An example function declaration */

__global__
void diagonal(int *arr,double *diag_array,int n,double coupling);

#endif // FOO_H_