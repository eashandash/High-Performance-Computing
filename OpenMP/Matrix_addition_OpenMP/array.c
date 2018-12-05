#include <stdio.h>
#include <omp.h>
#include <time.h>

int main()
{

  
   long double first[100000], second[100000];
   int c,d;
 
   
 /*
   for (c = 0; c < 200; c++)
   {
      for (d = 0; d < 200; d++)
      {
         first[c][d]=c;
         second[c][d]=d;
      }
   }
   */

 
   
   #ifdef _OPENMP
   float start=omp_get_wtime();
   #endif

   #pragma omp parallel for default(none) private(c) shared(first,second)
      for (c = 0; c < 100000; c++) 
      {
         first[c]=c+1;
         second[c]=(c+1)*2;
      }
         
    #ifdef _OPENMP 
    float end=omp_get_wtime();
    #endif     
      
    float result=end-start;

    printf("the time taken is %f\t",result);
    //printf("the time taken is 0.009937");


 
   return 0;
      }
