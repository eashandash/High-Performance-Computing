#include <stdio.h>
#include <omp.h>
#include <time.h>
#include<stdlib.h>

int main()
{

  
   int first[384][512], second[384][512], third[384][512];
   int c,d;
 int red[384][512], green[384][512], blue[384][512];

FILE* fpg,*fpr,*fpb;
fpg=fopen("green.txt","r");
fpr=fopen("red.txt","r");
fpb=fopen("blue.txt","r");

int i,j,k,l;

for(i=0;i<384;i++)
for(j=0;j<512;j++)
{
fscanf(fpg,"%d",&third[i][j]);
fscanf(fpr,"%d",&first[i][j]);
fscanf(fpg,"%d",&second[i][j]);
}


   //initialize the matrix
/* 
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
 
 for (c = 0; c < 400; c++)
   {
      for (d = 0; d < 200; d++)
      {
         //first[c][d]=c + 2*c + 2*d;
         third[c][d]=d+c;
      }
   }
 
 */

int threads[] = {1, 2, 4, 6, 8, 10, 12, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62};

 for(int m=0; m<20; m++)
{


 
   
   #ifdef _OPENMP
   float start=omp_get_wtime();
   #endif
   
   

   //#pragma omp parallel for default(none) shared(first,second,third,red,green,blue) private(i,j,k,l) 
   #pragma omp parallel for num_threads(threads[m]) shared(first,second,third,red,green,blue) private(i,j,k,l)         
         for(i = 0; i < 384; ++i)
            for(j = 0; j < 512; ++j)
	{
	if(i!=0 || i!=383 || j!=0 ||j!=511)
                        {	

	int sum=0;
                          for(k=i-1;k<=i+1;k++)
                            for(l=j-1;l<=j+1;l++)
	{
	sum=sum+first[k][l];

	}
	sum=sum/9;
	red[i][j]=sum;

	}


	}


for(i = 0; i < 384; ++i)
            for(j = 0; j < 512; ++j)
	{

	if(i!=0 || i!=383 || j!=0 ||j!=511)
                        {
	int sum1=0;
                          for(k=i-1;k<=i+1;k++)
                            for(l=j-1;l<=j+1;l++)
	{
	sum1=sum1+second[k][l];

	}
	sum1=sum1/9;
	blue[i][j]=sum1;

                        }


	}


for(i = 0; i < 384; ++i)
            for(j = 0; j < 512; ++j)
	{

	if(i!=0 || i!=383 || j!=0 ||j!=511)
                        {
	int sum2=0;
                          for(k=i-1;k<=i+1;k++)
                            for(l=j-1;l<=j+1;l++)
	{
	sum2=sum2+third[k][l];

	}
	sum2=sum2/9;
	green[i][j]=sum2;
	}	


	}





           

           
    #ifdef _OPENMP 
    float end=omp_get_wtime();
    #endif     
      
    float result=end-start;

    printf("the time taken for thread %d is %f\t \n",threads[m],result);
    
   /* for (c = 0; c < 5; c++)
   {
      for (d = 0; d < 5; d++)
      {
         printf(" %d ",mult[c][d]);
      }
      printf("\n");
      
   }
   
   */

  }
 

 
   return 0;
}
