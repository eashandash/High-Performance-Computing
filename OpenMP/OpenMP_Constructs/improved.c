#include <time.h>
#include <omp.h>
#include <stdio.h>



#define N 10
#define WIDTH 100
#define HEIGHT 100



int MAX;
int threads[12] = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 28, 32};

void Matrix(int M[WIDTH][HEIGHT])
{
	for (int i = 0; i < WIDTH; ++i)
	{
		for (int j = 0; j < HEIGHT; ++j)
		{
			M[i][j] = i+3*j;
		}
	}
}

void FindMax(int M[WIDTH][HEIGHT], int k)
{
	#pragma omp parallel for num_threads(threads[k]) shared(MAX)
	for (int i = 0; i < WIDTH; ++i)
	 {
			#pragma omp parallel for num_threads(threads[k]) shared(MAX)
			for (int j = 0; j < HEIGHT; ++j)
			{
				#pragma omp flush(MAX)
				if (M[i][j] > MAX)
				{
					#pragma omp critical
					if (M[i][j] > MAX)
					{
						MAX = M[i][j];
					}
				}
			}
	 } 	
}

int main(int argc, char const *argv[])
{
	int M1[WIDTH][HEIGHT], M2[WIDTH][HEIGHT], M3[WIDTH][HEIGHT], M4[WIDTH][HEIGHT], M5[WIDTH][HEIGHT], 
	M6[WIDTH][HEIGHT], M7[WIDTH][HEIGHT], M8[WIDTH][HEIGHT], M9[WIDTH][HEIGHT], M10[WIDTH][HEIGHT];
	
	//srand(time(0));

	Matrix(M1);
	Matrix(M2);
	Matrix(M3);
	Matrix(M4);
	Matrix(M5);
	Matrix(M6);
	Matrix(M7);
	Matrix(M8);
	Matrix(M9);
	Matrix(M10);

	float start, end, result;

	for (int k = 0; k < 12; ++k)
	{
		MAX = M1[0][0];

		start = omp_get_wtime();

		FindMax(M1, k);
		FindMax(M2, k);
		FindMax(M3, k);
		FindMax(M4, k);
		FindMax(M5, k);
		FindMax(M6, k);
		FindMax(M7, k);
		FindMax(M8, k);
		FindMax(M9, k);
		FindMax(M10, k);

		end = omp_get_wtime();
		result = end - start;
		printf("the time taken for thread %d is %f\t \n",threads[k],result);
	}
	//cout<<"Max element: "<<MAX<<endl;
	
	return 0;
}