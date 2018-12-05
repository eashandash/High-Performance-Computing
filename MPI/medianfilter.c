/******************************************************************************
* FILE: mpi_mm.c
* DESCRIPTION:  
*   MPI Matrix Multiply - C Version
*   In this code, the master task distributes a matrix multiply
*   operation to numtasks-1 worker tasks.
*   NOTE:  C and Fortran versions of this code differ because of the way
*   arrays are stored/passed.  C arrays are row-major order but Fortran
*   arrays are column-major order.
* AUTHOR: Blaise Barney. Adapted from Ros Leibensperger, Cornell Theory
*   Center. Converted to MPI: George L. Gusciora, MHPCC (1/95)
* LAST REVISED: 04/13/05
******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define N 200                   /* number of rows in matrix A */
#define M 200					/*number of columns in image*/	
#define MASTER 0                /* taskid of first task */
#define FROM_MASTER 1           /* setting a message type */
#define FROM_WORKER 2           /* setting a message type */

/*void swap(long double *p,long double *q) {
   long double t;
   
   t=*p; 
   *p=*q; 
   *q=t;
}

long double median(long double a[],int n) { 
   int i,j;

     for(i = 0;i < n-1;i++) {
        for(j = 0;j < n-i-1;j++) {
           if(a[j] > a[j+1])
              swap(&a[j],&a[j+1]);
        }
     }
   return(a[n/2]);
}*/

int main (int argc, char *argv[])
{
int	numtasks,              /* number of tasks in partition */
	taskid,                /* a task identifier */
	numworkers,            /* number of worker tasks */
	source,                /* task id of message source */
	dest,                  /* task id of message destination */
	mtype,                 /* message type */
	rows,                  /* rows of matrix A sent to each worker */
	averow, extra, offset, /* used to determine rows sent to each worker */
	i, j, k, rc;           /* misc */
long double	a[N][N],           /* matrix A IMAGE */
	b[N][N];           		   /* matrix B = RESULT OF MEAN FLTER */
MPI_Status status;

MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
if (numtasks < 2 ) {
  printf("Need at least two MPI tasks. Quitting...\n");
  MPI_Abort(MPI_COMM_WORLD, rc);
  exit(1);
  }
numworkers = numtasks-1;


/**************************** master task ************************************/
   if (taskid == MASTER)
   {
      printf("mpi_mm has started with %d tasks.\n",numtasks);
      printf("Initializing arrays...\n");
      for (i=0; i<N; i++)
         for (j=0; j<M; j++)
            a[i][j]= j;
            b[i][j]= j;

      /* Send matrix data to the worker tasks */
      averow = (N-2)/numworkers;
      extra = (N-2)%numworkers;
      printf("averow=%d, extra=%d \n",averow,extra);
      offset = 1;								//edges are ignored in mean filter algo
      mtype = FROM_MASTER;
      for (dest=1; dest<=numworkers; dest++)
      {
         rows = (dest <= extra) ? averow+1 : averow;   	
         printf("Sending %d rows to task %d offset=%d\n",rows,dest,offset);
         MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
         MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
         MPI_Send(&a[offset-1][0], (rows+2)*M, MPI_LONG_DOUBLE, dest, mtype,
                   MPI_COMM_WORLD);
         MPI_Send(&b[offset-1][0], (rows+2)*M, MPI_LONG_DOUBLE, dest, mtype, MPI_COMM_WORLD);
         offset = offset + rows;
      }
      printf("hellolefhailfj\n");
      /* Receive results from worker tasks */
      mtype = FROM_WORKER;
      for (i=1; i<=numworkers; i++)
      {
         source = i;
         MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
         MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
         MPI_Recv(&b[offset][0], rows*M, MPI_LONG_DOUBLE, source, mtype, 
                  MPI_COMM_WORLD, &status);
         printf("Received results from task %d\n",source);
      }

      /* Print results */
      printf("******************************************************\n");
      printf("Result Matrix:\n");
      for (i=0; i<N; i++)
      {
         printf("\n"); 
         for (j=0; j<N; j++) 
            printf("%6.2Lf   ", b[i][j]);
      }
      printf("\n******************************************************\n");
      printf ("Done.\n");
   }

/**************************** worker task ************************************/
   if (taskid > MASTER)
   {
        mtype = FROM_MASTER;
        MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&a, (rows+2)*M, MPI_LONG_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&b, (rows+2)*M, MPI_LONG_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);

        long double c[9],temp;
        int p,q;
        //printf("taskid=%d\n", taskid);
        for (i=1; i<rows+1; i++)
        {	
		        for (j=1; j<M-1; j++)
	           {   
                 c[0]=a[i-1][j-1] ; c[1]=a[i][j-1] ; c[2]=a[i+1][j-1] ; c[3]=a[i+1][j] ; c[4]=a[i+1][j+1] ; c[5]=a[i][j+1] ; c[6]=a[i-1][j+1] ; c[7]=a[i-1][j] ;c[8]=b[i][j];
	                  
                 for(p = 0;p < 8;p++) {
                   for(q = 0;q < 9-i-1;q++) {
                      if(c[j] > c[j+1])
                        temp=c[j];
                        c[j+1]=temp;
                        c[j]=c[j+1]; 
                      }
                    }
                 
                 b[i][j] = c[4];
	           }
	      }	
        mtype = FROM_WORKER;
        MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        MPI_Send(&b[1], rows*M, MPI_LONG_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
  }
   MPI_Finalize();
}
