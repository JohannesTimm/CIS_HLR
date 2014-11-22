/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>

#include "partdiff-seq.h"


#ifndef MODE
#define MODE 0  // 0 = sequentiell, 1 = mit Posix Threads
#endif
#if MODE==1
#include <pthread.h>
#endif
#define NUM_THREADS 12

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};


typedef struct t_thread_data
{
   long  thread_id;
   struct calculation_arguments *arguments;
   struct calculation_results *results;
   struct options *options;
}t_thread_data;

struct t_thread_data thread_data_array[NUM_THREADS] ;

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes)\n", size);
		/* exit program */
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
/*  local variables for loops   */
  uint64_t g, i, j;

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}
# if (MODE==1)
  void *calculate_wrapper(void * params) {
    long id;
    t_thread_data *p = (t_thread_data *)params;
    struct calculation_arguments *arguments = p->arguments;
    
    // Das geht nicht, warum??
    // calculate (p->arguments, p->results, p->options);
    return NULL;
  }
#endif
typedef struct  //give parameters for every thread
{
    struct calculation_arguments const* arguments;
    struct options const* options;
    
    int thread_id;
    int m1; //the matrixes
    int m2;
    int i_start;
    int i_end;
    double* residuum;
    double* maxresiduum;
    int term_iteration;
} thread_args;

static pthread_mutex_t mutex_residuum= PTHREAD_MUTEX_INITIALIZER;
/* ************************************************************************ */
/*threaded_calc : */
/* ************************************************************************ */
static
void*
threaded_calc(void* params)
{
	
	 int i, j;
	double star;
	thread_args* args= (thread_args*) params;
	double** Matrix_Out = args->arguments->Matrix[args->m1];
	double** Matrix_In  = args->arguments->Matrix[args->m2];
	int i_start= args->i_start;
	int i_end=args ->i_end;
	if (args->options->inf_func == FUNC_FPISIN)
	{
	double* fpisin_i;
	double pih;
	double fpisin;
	fpisin_i = (double*)malloc(args->arguments->N * sizeof(double));
	pih = PI * args->arguments->h;
	fpisin = 0.25 * TWO_PI_SQUARE * args->arguments->h * args->arguments->h;
	for (i = i_start; i < i_end; i++)			
				{ 	
					fpisin_i[i] = 0.0;
					fpisin_i[i]= fpisin * sin(pih * (double)i);
				}
	}
	/* over all rows */
	for (i = i_start; i < i_end; i++)		
 	{
		/* over all columns */
        for (j = 1; j < args->arguments->N; j++)
        {
		star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

					if (args->options->inf_func == FUNC_FPISIN)
					{
						star += fpisin_i[i] * sin(pih * (double)j);
					}
		if (args->options->termination == TERM_PREC || args->term_iteration == 1)
            {
                pthread_mutex_lock(&mutex_residuum); *args->residuum = Matrix_In[i][j] - star;
                *args->residuum = (*args->residuum < 0) ? -*args->residuum : *args->residuum;
                *args->maxresiduum = (*args->residuum < *args->maxresiduum) ? *args->maxresiduum : *args->residuum;
                pthread_mutex_unlock(&mutex_residuum);
            }
		Matrix_Out[i][j] = star;
		}
	}
}










/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options)
{
	int i, j;                                   /* local variables for loops  */
	int m1, m2;                                 /* used as indices for old and new matrices       */
	//double star;                                /* four times center value minus 4 neigh.b values */
        double res = 0;
        double maxres = 0;
	double* residuum = &res;                            /* residuum of current iteration                  */
	double* maxresiduum = &maxres;                         /* maximum residuum value of a slave in iteration */


	int nthreads= options->number; //add here the number from the options

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	pthread_t threads[nthreads]; //build our threads
//**
//Aufteilungsalgo here!!!!!
	int width;
	int rest;

	width = (int)((N-1)/nthreads);
	rest = (N-1) % nthreads;
	
	int istart[nthreads];
	int iend[nthreads];
	
	for (i=0; i<nthread; i++)
	{
		if (i < rest)
		{
			istart[i] = (width + 1) * i + 1;
			iend[i] = istart[i] + width + 1;
		}
		else
		{
			istart[i] = width * i + rest + 1;
			iend[i] = istart[i] + width;
		}
	}
//+++
//+++
//+++
//+++
//+++

	while (term_iteration > 0)
	{
		maxresiduum = 0;
		
		//creating threads here
		for(i = 0; i < nthreads; ++i)
        {	
			int rc; //test if sucessful
            thread_args* args = malloc(sizeof(thread_args));
            args->arguments = arguments;
            args->options = options;
            args->thread_id = i;
            args->m1 = m1;
            args->m2 = m2;
            args->i_start=istart[i];
            args->i_end=iend[i];
            args->residuum = residuum;
            args->maxresiduum = maxresiduum;
            args->term_iteration = term_iteration;
            rc=pthread_create(&threads[i], NULL, threaded_calc, (void *)args);
			if (rc)
			{
	         	printf("ERROR; return code from pthread_create() is %d\n", rc);
	        	exit(-1);
			}
		 }
		  //joining them
		for(i = 0; i < nthreads; ++i)
        {
			int rc;			            
			rc=pthread_join(threads[i], NULL);
	
    		  if (rc) 
				{
    		     printf("ERROR; return code from pthread_join() is %d\n", rc);
         		 exit(-1);
         		}
        }
					
    
        


//		double** Matrix_Out = arguments->Matrix[m1];
//		double** Matrix_In  = arguments->Matrix[m2];
//
//		maxresiduum = 0;
//
//		/* over all rows */
//		for (i = 1; i < N; i++)
//		{
//			double fpisin_i = 0.0;
//
//			if (options->inf_func == FUNC_FPISIN)
//			{
//				fpisin_i = fpisin * sin(pih * (double)i);
//			}
//
//			/* over all columns */
//			for (j = 1; j < N; j++)
//			{
//				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);
//
//				if (options->inf_func == FUNC_FPISIN)
//				{
//					star += fpisin_i * sin(pih * (double)j);
//				}
//
//				if (options->termination == TERM_PREC || term_iteration == 1)
//				{
//					residuum = Matrix_In[i][j] - star;
//					residuum = (residuum < 0) ? -residuum : residuum;
//					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
//				}
//
//				Matrix_Out[i][j] = star;
//			}
//		}

		results->stat_iteration++;
		results->stat_precision = *maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation, depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (*maxresiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauss-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion DisplayMatrix:                               **/
/**                                                                        **/
/** Die Funktion DisplayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;
# if (MODE==1)
  long t;  // thread number
  int rc;  // return value
  struct thread_data *data;
  pthread_t threads[NUM_THREADS];
#endif

/* ************************* */
/* get parameters */
/* ************************* */
  AskParams(&options, argc, argv);
  #if (MODE==0)
    // 0 = sequentieller Code
	  printf("Sequentieller Code\n");
    // ...

  #elif (MODE == 1)
    // 1 = Aufteilung mit Posix-Threads
    printf("Aufteilung auf %d Posix Threads\n", (int) options.number);
    if (options.number > NUM_THREADS) {
      printf("ERROR; Mehr als %d Threads geht nicht!\n", NUM_THREADS);
      exit (-1);
    }
    for (t = 0; t < (long) options.number; t++) {
      rc = pthread_create(&threads[t], NULL, calculate_wrapper, (void *)t);
      if (rc){
        printf("ERROR; return code from pthread_create() is %d\n", rc);
        exit(-1);
      }
      
    }
  pthread_exit(NULL);
  #else
	  printf("Ungültiger Mode = %d\n", MODE); // Nur zum Testen
  #endif

/* ******************************************* */
/*  get and initialize variables and matrices  */
/* ******************************************* */
	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

/*  start timer         */
	gettimeofday(&start_time, NULL);
/*  solve the equation  */
	calculate(&arguments, &results, &options);
/*  stop timer          */
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	DisplayMatrix(&arguments, &results, &options);

/*  free memory     */
	freeMatrices(&arguments);

	return 0;
}
