#include <stdio.h>
#include <omp.h>
#include <time.h>

int main()
{

  
   int first[200][400], second[400][200], mult[200][200];
   int c,d;
 
   //initialize the matrix
 
   for (c = 0; c < 200; c++)
   {
      for (d = 0; d < 400; d++)
      {
         first[c][d]=c + 2*c + 2*d;
         //second[c][d]=d;
      }
   }
 
   for (c = 0; c < 400; c++)
   {
      for (d = 0; d < 200; d++)
      {
         //first[c][d]=c + 2*c + 2*d;
         second[c][d]=d;
      }
   }
 
 
 
 
   
   #ifdef _OPENMP
   float start=omp_get_wtime();
   #endif
   
   int i,j,k;

   #pragma omp parallel for default(none) private(i,j,k) shared(mult,first,second)            
         for(j = 0; j < 200; ++j)
            for(i = 0; i < 200; ++i)
            for(k = 0; k < 400; ++k)
            {
                mult[i][j] += first[i][k] * second[k][j];
            }

           
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
