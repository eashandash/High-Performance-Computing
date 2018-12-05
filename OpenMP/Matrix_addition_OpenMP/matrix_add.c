#include <stdio.h>
#include <omp.h>
#include <time.h>

int main()
{

  
   long double first[200][200], second[200][200], sum[200][200];
   int c,d;
 
   
 
   for (c = 0; c < 200; c++)
   {
      for (d = 0; d < 200; d++)
      {
         first[c][d]=c;
         second[c][d]=d;
      }
   }
 
   
   #ifdef _OPENMP
   float start=omp_get_wtime();
   #endif

   #pragma omp parallel for default(none) private(c,d) shared(sum,first,second)
      for (c = 0; c < 200; c++) 
         for (d = 0 ; d < 200; d++) 
            sum[c][d] = first[c][d] + second[c][d];

    #ifdef _OPENMP 
    float end=omp_get_wtime();
    #endif     
      
    float result=end-start;

    printf("the time taken is %f\t",result);

 
   return 0;
}