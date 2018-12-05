#include <stdio.h>
#include <omp.h>
#include <time.h>

#define N 100

#define B 8

int min(int a,int b)
{
if (a<=b)
return a;

else
return b;

}




int main()
{

  
  
   int y[N][N], z[N][N], x[N][N];
   int c,d;
 
   //initialize the matrix
 
   for (c = 0; c < N; c++)
   {
      for (d = 0; d < N; d++)
      {
         y[c][d]=c + 2*c + 2*d;
         //second[c][d]=d;
      }
   }
 
   for (c = 0; c < N; c++)
   {
      for (d = 0; d < N; d++)
      {
         //first[c][d]=c + 2*c + 2*d;
         z[c][d]=d;
      }
   }
 
 
 
 
   
   #ifdef _OPENMP
   float start=omp_get_wtime();
   #endif
   
   int i,j,k;
   int r;
   int kk,jj;

   #pragma omp parallel for default(none) private(i,j,k,kk,jj) shared(x,y,z,r)          
         for (jj = 0; jj < N; jj = jj+B)
          for (kk = 0; kk < N; kk = kk+B)
           for (i = 0; i < N; i = i+1)
            for (j = jj; j < min(jj+B,N); j = j+1)
             {r = 0;
             for (k = kk; k < min(kk+B,N); k = k+1) {
             r = r + y[i][k]*z[k][j];};
             x[i][j] = x[i][j] + r;
             };

           
    #ifdef _OPENMP 
    float end=omp_get_wtime();
    #endif     
      
    float result=end-start;

    printf("the time taken is %f\t \n",result);
   


   /* for (c = 0; c < 5; c++)
   {
      for (d = 0; d < 5; d++)
      {
         printf(" %d ",mult[c][d]);
      }
      printf("\n");
      
   }
   
   */
 

 
   return 0;
}
