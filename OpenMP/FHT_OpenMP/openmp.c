
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>




#include <math.h>
#if (_OPENMP >= 200203)
    #include <omp.h>
#endif



#ifndef FFT_UNIT_STRIDE
    #define FFT_UNIT_STRIDE 0
#endif
#ifndef FHT_UNIT_STRIDE
    #define FHT_UNIT_STRIDE 0
#endif



#define PI 3.14
#define L1_CACHE_BYTES 2

static void fht_dif_iter(double*, unsigned long);
static void fht_dif_iter_seq(double*, unsigned long);
static void fht_dif_rec(double*, unsigned long, int);

static void fht_dit_iter(double*, unsigned long);
static void fht_dit_iter_seq(double*, unsigned long);
static void fht_dit_rec(double*, unsigned long, int);



/***********************************************************
* void fht_dif(double *x, unsigned long n)
*
* Purpose:
*   Computes the discrete Hartley transform of a real sequence
*   x[0..n-1], using a radix 2 decimation-in-frequency FHT.
*   n must be a power of 2.  Entering this function, x[] is
*   in normal order.  On return, x[] contains the Hartley
*   transform, stored in bit-reversed order.
*
* Arguments:
*   double *x - array of n doubles, representing n real numbers
*   unsigned long n - dimension of x, must be a power of 2
*
* Return:
*   none
************************************************************/

void fht_dif(double *x, unsigned long n)
{
    fht_dif_rec(x, n, 1);
    return;
}   


/***********************************************************
* void fht_dit(double *x, unsigned long n)
*
* Purpose:
*   Computes the discrete Hartley transform of a real sequence
*   x[0..n-1], using a radix 2 decimation-in-time FHT.
*   n must be a power of 2.  Entering this function, x[]
*   must be in bit-reversed order.  On return, x[] contains
*   the Hartley transform, returned to normal order.
*
* Arguments:
*   double *x - array of n doubles, representing n real numbers
*   unsigned long n - dimension of x, must be a power of 2
*
* Return:
*   none
************************************************************/

/***********************************************************
* static void fht_dif_iter(double *x, unsigned long n)
*
* Purpose:
*   Computes the discrete Hartley transform of a real sequence
*   x[0..n-1], using an iterative radix 2 decimation-in-frequency
*   FHT.  n must be a power of 2.  Entering this function, x[] is
*   in normal order.  On return, x[] contains the Hartley
*   transform, stored in bit-reversed order.
*
* Arguments:
*   double *x - array of n doubles, representing n real numbers
*   unsigned long n - dimension of x, must be a power of 2
*
* Return:
*   none
************************************************************/
int main(int argc, char const *argv[])
{
	
	//printf("ED\n");
	
	double a[16] = {1,2,3,4,5,6,7,8,9,12,13,14,15,16,17,18};

	unsigned long n = 16;
	int i;






	//fht_dif_iter(a, n);
	printf("\n Welcome to Discrete Hartley Transform   \n");
	


	printf("\n The 1-D signal under consideration is  :  \n ");
	for(i=0;i<16;i++)
{
   printf(" %lf ",a[i]);
}

printf("\n");

printf("\n The Discrete Hartley frequency components are : \n");



	
	// fht_dif_iter_seq(a, n);
	
	//fht_dif_rec(a, n, 10);
	
	// fht_dit_iter(a, n);
	
	// fht_dit_iter_seq(a, n);


	
	  
   #ifdef _OPENMP
   float start=omp_get_wtime();
   #endif


	fht_dit_rec(a, n,10);


	  #ifdef _OPENMP 
    float end=omp_get_wtime();
    #endif     
      
    float result=end-start;

    printf("the time taken is %f\t",result);

	 

	return 0;
}

//function definations


static void fht_dif_iter(double *x, unsigned long n)
{
    unsigned long m;

    for (m = n; m > 1; m >>= 1) {
        double a, b, c, s, t;
        unsigned long i, j, k, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (i = 0; i < n; i += m) {
            double *xp;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k) {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;

                // printf("j : %lu \t %f\n",j, xp[j]);
                // printf("k : %lu \t %f\n",k, xp[k]);


                                 
            }
        }
        c = 1.0;
        s = 0.0;
        for (j = 1, k = mh - 1; j < mq; ++j, --k) {
            double tmp;
            double *xj, *xk;
            xj = x + j + mh;
            xk = x + k + mh;
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
            for (i = 0; i < n; i += m) {
                double u, v;
                u = xj[i];
                v = xk[i];
                xj[i] = u * c + v * s;
                xk[i] = u * s - v * c;


                
            }
        }
    }

   
    return;
}  





/***********************************************************
* static void fht_dit_rec(double *x, unsigned long n, int nbranch)
*
* Purpose:
*   Computes the discrete Hartley transform of a real sequence
*   x[0..n-1], using a radix 2 decimation-in-time FHT.  If
*   the computation is small enough to fit in cache, it is done
*   iteratively.  Otherwise, it is done recursively until the
*   recursion descends to cache-sized chunks.
*
*   n must be a power of 2.  Entering this function, x[] must
*   be in bit-reversed order.  On return, x[] contains the
*   Hartley transform, returned to normal order.
*
*   To support OpenMP parallelism, nbranch keeps track of the
*   number of active transforms at a given recursion level. On
*   the first call to this function, nbranch should be 1.  It
*   is then doubled for each recursion.
*
* Arguments:
*   double *x - array of n doubles, representing n real numbers
*   unsigned long n - dimension of x, must be a power of 2
*   int nbranch - number of transforms at this recursion level
*
* Return:
*   none
************************************************************/

static void fht_dit_rec(double *x, unsigned long n, int nbranch)
{
    double a, b, c, s, t;
    unsigned long j, jmax, k, nh, nq;

    if (n == 1)
        return;
    if (n <= (unsigned long)(L1_CACHE_BYTES / sizeof(double))) {
        if (FHT_UNIT_STRIDE)
            fht_dit_iter_seq(x, n);
        else
            fht_dit_iter(x, n);
        return;
    }
    nh = n >> 1;
    nq = nh >> 1;
    nbranch <<= 1;

     


    #if (_OPENMP >= 200203)
    #pragma omp parallel sections if (nbranch <= omp_get_max_threads()) num_threads(2)
    #endif
    {
        #if (_OPENMP >= 200203)
        #pragma omp section
        #endif
        fht_dit_rec(x, nh, nbranch);
        #if (_OPENMP >= 200203)
        #pragma omp section
        #endif
        fht_dit_rec(x + nh, nh, nbranch);
    }




    t = PI / (double)nh;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    jmax = nq + nh;
    c = 1.0;
    s = 0.0;
    for (j = nh + 1, k = n - 1; j < jmax; ++j, --k) {
        double tmp, u, v;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
        u = x[j];
        v = x[k];
        x[j] = u * c + v * s;
        x[k] = u * s - v * c;

         
    }
    for (j = 0, k = nh; j < nh; ++j, ++k) {
        double u, v;
        u = x[j];
        v = x[k];
        x[j] = u + v;
        x[k] = u - v;

        // printf("j : %lu \t %f\n",j, x[j]);
        // printf("k : %lu \t %f\n",k, x[k]);
        printf("%f\n,",x[j]- x[k]);
    }


 
  

    return;
}   




/***********************************************************
* static void fht_dif_rec(double *x, unsigned long n, int nbranch)
*
* Purpose:
*   Computes the discrete Hartley transform of a real sequence
*   x[0..n-1], using a radix 2 decimation-in-frequency FHT.  If
*   the computation is small enough to fit in cache, it is done
*   iteratively.  Otherwise, it is done recursively until the
*   recursion descends to cache-sized chunks.
*
*   n must be a power of 2.  Entering this function, x[] must
*   be in normal order.  On return, x[] contains the Hartley
*   transform, in bit-reversed order.
*
*   To support OpenMP parallelism, nbranch keeps track of the
*   number of active transforms at a given recursion level. On
*   the first call to this function, nbranch should be 1.  It
*   is then doubled for each recursion.
*
* Arguments:
*   double *x - array of n doubles, representing n real numbers
*   unsigned long n - dimension of x, must be a power of 2
*   int nbranch - number of transforms at this recursion level
*
* Return:
*   none
************************************************************/

static void fht_dif_rec(double *x, unsigned long n, int nbranch)
{
    double a, b, c, s, t;
    unsigned long j, jmax, k, nh, nq;

    if (n == 1)
        return;
    if (n <= (unsigned long)(L1_CACHE_BYTES / sizeof(double))) {
        if (FHT_UNIT_STRIDE)
            fht_dif_iter_seq(x, n);
        else
            fht_dif_iter(x, n);
        return;
    }
    nh = n >> 1;
    nq = nh >> 1;
    t = PI / (double)nh;
    a = sin(0.5 * t);
    a *= 2.0 * a;
    b = sin(t);
    for (j = 0, k = nh; j < nh; ++j, ++k) {
        double u, v;
        u = x[j];
        v = x[k];
        x[j] = u + v;
        x[k] = u - v;
    }
    c = 1.0;
    s = 0.0;
    jmax = nq + nh;
    for (j = nh + 1, k = n - 1; j < jmax; ++j, --k) {
        double u, v, tmp;
        tmp = c;
        c -= a * c + b * s;
        s -= a * s - b * tmp;
        u = x[j];
        v = x[k];
        x[j] = u * c + v * s;
        x[k] = u * s - v * c;

       

    }
    nbranch <<= 1;
    #if (_OPENMP >= 200203)
    #pragma omp parallel sections if (nbranch <= omp_get_max_threads()) num_threads(2)
    #endif
    {
        #if (_OPENMP >= 200203)
        #pragma omp section
        #endif
        fht_dif_rec(x, nh, nbranch);
        #if (_OPENMP >= 200203)
        #pragma omp section
        #endif
        fht_dif_rec(x + nh, nh, nbranch);
    }
    return;
}  














/***********************************************************
* static void fht_dif_iter_seq(double *x, unsigned long n)
*
* Purpose:
*   Computes the discrete Hartley transform of a real sequence
*   x[0..n-1], using an iterative radix 2 decimation-in-frequency
*   FHT.  n must be a power of 2.  Entering this function, x[] is
*   in normal order.  On return, x[] contains the Hartley
*   transform, stored in bit-reversed order.
*
*   The two inner loops of the FHT computation are ordered
*   to favor sequential memory access at the expense of
*   redundant trig computations.  See J. Arndt, "Algorithms
*   for Programmers," online at http://www.jjj.de/fxt/.
*
* Arguments:
*   double *x - array of n doubles, representing n real numbers
*   unsigned long n - dimension of x, must be a power of 2
*
* Return:
*   none
************************************************************/

static void fht_dif_iter_seq(double *x, unsigned long n)
{
    unsigned long m;

    for (m = n; m > 1; m >>= 1) {
        double a, b, t;
        unsigned long i, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (i = 0; i < n; i += m) {
            double c, s;
            double *xp;
            unsigned long j, k;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k) {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
            xp += mh;
            c = 1.0;
            s = 0.0;
            for (j = 1, k = mh - 1; j < mq; ++j, --k) {
                double u, v, tmp;
                tmp = c;
                c -= a * c + b * s;
                s -= a * s - b * tmp;
                u = xp[j];
                v = xp[k];
                xp[j] = u * c + v * s;
                xp[k] = u * s - v * c;
            }
        }
    }
    return;
}   




/***********************************************************
* static void fht_dit_iter(double *x, unsigned long n)
*
* Purpose:
*   Computes the discrete Hartley transform of a real sequence
*   x[0..n-1], using an iterative radix 2 decimation-in-time
*   FHT.  n must be a power of 2.  Entering this function, x[]
*   must be in bit-reversed order.  On return, x[] contains the
*   Hartley transform, returned to normal order.
*
* Arguments:
*   double *x - array of n doubles, representing n real numbers
*   unsigned long n - dimension of x, must be a power of 2
*
* Return:
*   none
************************************************************/

static void fht_dit_iter(double *x, unsigned long n)
{
    unsigned long m;

    for (m = 2; m <= n; m <<= 1) {
        double a, b, c, s, t;
        unsigned long i, j, k, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        c = 1.0;
        s = 0.0;
        for (j = 1, k = mh - 1; j < mq; ++j, --k) {
            double tmp;
            double *xj, *xk;
            xj = x + j + mh;
            xk = x + k + mh;
            tmp = c;
            c -= a * c + b * s;
            s -= a * s - b * tmp;
            for (i = 0; i < n; i += m) {
                double u, v;
                u = xj[i];
                v = xk[i];
                xj[i] = u * c + v * s;
                xk[i] = u * s - v * c;
            }
        }
        for (i = 0; i < n; i += m) {
            double *xp;
            xp = x + i;
            for (j = 0, k = mh; j < mh; ++j, ++k) {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
        }
    }
    return;
}   


/***********************************************************
* static void fht_dit_iter_seq(double *x, unsigned long n)
*
* Purpose:
*   Computes the discrete Hartley transform of a real sequence
*   x[0..n-1], using an iterative radix 2 decimation-in-time
*   FHT.  n must be a power of 2.  Entering this function, x[]
*   must be in bit-reversed order.  On return, x[] contains the
*   Hartley transform, returned to normal order.
*
*   The two inner loops of the FHT computation are ordered
*   to favor sequential memory access at the expense of
*   redundant trig computations.  See J. Arndt, "Algorithms
*   for Programmers," online at http://www.jjj.de/fxt/.
*
* Arguments:
*   double *x - array of n doubles, representing n real numbers
*   unsigned long n - dimension of x, must be a power of 2
*
* Return:
*   none
************************************************************/

static void fht_dit_iter_seq(double *x, unsigned long n)
{
    unsigned long m;

    for (m = 2; m <= n; m <<= 1) {
        double a, b, t;
        unsigned long i, mh, mq;
        mh = m >> 1;
        mq = mh >> 1;
        t = PI / (double)mh;
        a = sin(0.5 * t);
        a *= 2.0 * a;
        b = sin(t);
        for (i = 0; i < n; i += m) {
            double c, s;
            double *xp;
            unsigned long j, k;
            xp = x + i + mh;
            c = 1.0;
            s = 0.0;
            for (j = 1, k = mh - 1; j < mq; ++j, --k) {
                double tmp, u, v;
                tmp = c;
                c -= a * c + b * s;
                s -= a * s - b * tmp;
                u = xp[j];
                v = xp[k];
                xp[j] = u * c + v * s;
                xp[k] = u * s - v * c;
            }
            xp -= mh;
            for (j = 0, k = mh; j < mh; ++j, ++k) {
                double u, v;
                u = xp[j];
                v = xp[k];
                xp[j] = u + v;
                xp[k] = u - v;
            }
        }
    }
    return;
}   

