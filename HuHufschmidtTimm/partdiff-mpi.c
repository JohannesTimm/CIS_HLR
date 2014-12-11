 /****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-mpi.c                                              **/
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
#include <mpi.h>

#include "partdiff-mpi.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
	int N_local;
	int rank;
	int size;
	int from;
	int to;
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

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
	int rest;
	uint64_t const N = arguments->N;
	
	MPI_Comm_rank (MPI_COMM_WORLD, &arguments -> rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &arguments -> size);
	
	int const size = arguments -> size;
	int const rank = arguments -> rank;	
	
	arguments->N = (options->interlines * 8) + 9 - 1;
	rest = (N + 1 -2) % size;
		
	if(rank < rest)
	{
		arguments -> N_local = (N + 1 -2) / size + 1;
		arguments -> from = rank * ((N - 1)/size + 1) + 1;
		arguments -> to = arguments -> from + (N-1)/size;
	}
	else 
	{
		arguments -> N_local = (N + 1 -2) / size;
		arguments -> from = rank * ((N - 1)/size) + rest + 1;
		arguments -> to = arguments -> from + (N-1)/size - 1;
	}
	
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
	int const N_local = arguments -> N_local;

	arguments->M = allocateMemory(arguments->num_matrices * (N_local + 2) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N_local + 2) * sizeof(double*));

		for (j = 0; j < N_local + 2 ; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N_local + 2) * (N + 1)) + (j * (N + 1));
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
	uint64_t g, i, j;                                /*  local variables for loops   */

	uint64_t const N = arguments -> N;
	int const N_local =arguments -> N_local;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;
	int from = arguments -> from;
	int rank = arguments -> rank;
	int size = arguments -> size;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i < N_local + 2; i++)
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
			for (i = 0; i < N_local + 2; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i) - (from - 1) * h;
				Matrix[g][i][N] = h * i + (from - 1) * h;
			}
			if (rank == 0)	
			{	
				for (i = 0; i <= N; i++)
				{
					Matrix[g][0][i] = 1.0 - (h * i);
				}
			}
			if (rank + 1 == size)
			{
				for (i = 0; i <= N; i++)
				{
					Matrix[g][N][i] = h * i;
				}
			}

			Matrix[g][N_local + 1][0] = 0.0;
			Matrix[g][0][N] = 0.0;
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
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration                  */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

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

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = 0;

		/* over all rows */
		for (i = 1; i < N; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation, depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxresiduum < options->term_precision)
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


///*************************************************************************** */
///* calculateMPI: calculates using MPI communication
///* ************************************************************************** */
//static
//void
//calculate (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options)
//{
	//int i, j;                                   /* local variables for loops  */
	//int m1, m2;                                 /* used as indices for old and new matrices       */
	//double star;                                /* four times center value minus 4 neigh.b values */
	//double residuum;                            /* residuum of current iteration                  */
	//double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	//int const N = arguments->N;
	//double const h = arguments->h;

	//double pih = 0.0;
	//double fpisin = 0.0;

	//int term_iteration = options->term_iteration;

	///* initialize m1 and m2 depending on algorithm */
	//if (options->method == METH_JACOBI)
	//{
		//m1 = 0;
		//m2 = 1;
	//}
	//else
	//{
		//m1 = 0;
		//m2 = 0;
	//}

	//if (options->inf_func == FUNC_FPISIN)
	//{
		//pih = PI * h;
		//fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	//}

	//while (term_iteration > 0)
	//{
		//double** Matrix_Out = arguments->Matrix[m1];
		//double** Matrix_In  = arguments->Matrix[m2];

		//maxresiduum = 0;
		
		////openmp here
		////openmp reduce max maxresiduum
		///* over all rows */
		//for (i = 1; i < N_local; i++)
		//{
			//double fpisin_i = 0.0;

			//if (options->inf_func == FUNC_FPISIN)
			//{
				//fpisin_i = fpisin * sin(pih * (double)i);
			//}
			
			//if (i=1) //for first line excahnge halos (onlz the thread actuallz computing this has toit in the jacobi scheme)
			//{ 
				//if (rank>0) //exclude top rank
				//{
					//MPI_SENDRECV()
				//} 
			//}
			
			//if (i=N_local-1) // exchange last halo lines
			//{	
				//if (rank<max_rank) // exclude last rank
				//{
					//MPI_SENDRECV()
				//}
			//}
			 
			///* over all columns */
			//for (j = 1; j < N; j++)
			//{
				//star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				//if (options->inf_func == FUNC_FPISIN)
				//{
					//star += fpisin_i * sin(pih * (double)j);
				//}

				//if (options->termination == TERM_PREC || term_iteration == 1)
				//{
					//residuum = Matrix_In[i][j] - star;
					//residuum = (residuum < 0) ? -residuum : residuum;
					//maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				//}

				//Matrix_Out[i][j] = star;
			//}
		//}
		////end openmp
		//localmaxresiduum // got form openmp

		////this needs to be thread/mpi safe!!
		//results->stat_iteration++;
		//results->stat_precision = maxresiduum;

		///* exchange m1 and m2 */
		//i = m1;
		//m1 = m2;
		//m2 = i;

		///* check for stopping calculation, depending on termination method */
		//if (options->termination == TERM_PREC)
		//{
			//MPIALLREDUCE (Max GlobalMax, localMax)
			//if (rank=0)
			//{
				//if (maxresiduum < options->term_precision)
				//{
				//term_iteration = 0;
				//}
				//else
				//{ 
				//term_iteration =1;
				//}
			//}
			//MPI_BCAST(from rank 0 to everzbodz, term_iteration,1,MPI_INT)	
		//}
		//else if (options->termination == TERM_ITER) // works without special syncronisation
		//{
			//term_iteration--;
		//}
	//}

	//results->m = m2;
//}
		
		
///* ************************************************************************* */
///* sub_iHALOExcahnage: controls the functions taht send halo lines         */
///* ************************************************************************ */
//static
//void
//sub_iHaloExchance(Piece,rank,numProcess,root,COMM,requestUp, requestDown)
	//{
		
		//double** Piece;
		//int root, COMM, numProcess, rank;
		
		//MPI_Request* requestUp, requestDown; //request[2]
		//MPI_Status   status;
		
		////!UP: Values, that goes from rank to rank-1 (If Root is on top, it goes UP)
		////! Tag = 100 + source
		////call sub_iHaloSendUp(Piece,rank,numProcess,root,MPI_COMM_WORLD,requestUp(1))
		//if (rank>root)
		//{
		 //if (MPI_WAIT(requestUp[1], status))!=MPI_SUCSESS)
		  //{printf("ERROR WAIT");}
		//} 
		
		////call sub_iHaloRecvUp(Piece,rank,numProcess,root,MPI_COMM_WORLD,requestUp(2))
		//if (rank<numProcess-1)
		//{
		//if (MPI_WAIT(requestUp[2], status))!=MPI_SUCSESS)
		  //{printf("ERROR WAIT");}
		 //}
		////!-------------------------------------------------------------------------------
		////!DOWN: Values, that goes from rank to rank+1 (If Root is on top, it goes DOWN) 
		////! Tag = 200 + source
		////call sub_iHaloSendDown(Piece,rank,numProcess,root,MPI_COMM_WORLD,requestDown(1))
		//if (rank<numProcess-1) 
		//{
		//if (MPI_WAIT(requestDown[1], status))!=MPI_SUCSESS)
		  //{printf("ERROR WAIT");}
		//}
		////call sub_iHaloRecvDown(Piece,rank,numProcess,root,MPI_COMM_WORLD,requestDown(2))
		//if (rank>root) 
		//{
		//if (MPI_WAIT(requestDown[2], status))!=MPI_SUCSESS)
		  //{printf("ERROR WAIT");}
		//}
	//}		
	
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

///****************************************************************************/
///** Beschreibung der Funktion DisplayMatrix:                               **/
///**                                                                        **/
///** Die Funktion DisplayMatrix gibt eine Matrix                            **/
///** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
///**                                                                        **/
///** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
///** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
///** sieben Zwischenzeilen ausgegeben.                                      **/
///****************************************************************************/
//static
//void
//DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
//{
	//int x, y;

	//double** Matrix = arguments->Matrix[results->m];

	//int const interlines = options->interlines;

	//printf("Matrix:\n");

	//for (y = 0; y < 9; y++)
	//{
		//for (x = 0; x < 9; x++)
		//{
			//printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		//}

		//printf ("\n");
	//}

	//fflush (stdout);
//}
/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are halo lines for other processes. It is responsible for (global) lines 5-6.
 */
static
void
//DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
  int const elements = 8 * options->interlines + 9;

  int x, y;
  double** Matrix = arguments->Matrix[results->m];
  MPI_Status status;
  
  int from = arguments -> from;
  int to = arguments -> to;
  int const size = arguments -> size;
  int const rank = arguments -> rank;
  /* first line belongs to rank 0 */
  if (rank == 0)
    from--;

  /* last line belongs to rank size - 1 */
  if (rank + 1 == size)
    to++;

  if (rank == 0)
    printf("Matrix:\n");

  for (y = 0; y < 9; y++)
  {
    int line = y * (options->interlines + 1);

    if (rank == 0)
    {
      /* check whether this line belongs to rank 0 */
      if (line < from || line > to)
      {
        /* use the tag to receive the lines in the correct order
         * the line is stored in Matrix[0], because we do not need it anymore */
        MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        /* if the line belongs to this process, send it to rank 0
         * (line - from + 1) is used to calculate the correct local address */
        MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
      }
    }

    if (rank == 0)
    {
      for (x = 0; x < 9; x++)
      {
        int col = x * (options->interlines + 1);

        if (line >= from && line <= to)
        {
          /* this line belongs to rank 0 */
          printf("%7.4f", Matrix[line][col]);
        }
        else
        {
          /* this line belongs to another rank and was received above */
          printf("%7.4f", Matrix[0][col]);
        }
      }

      printf("\n");
    }
  }

  fflush(stdout);
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
	int rc;
	
	/* mpi starts*/
	rc= MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) 
	{
		printf ("Error starting MPI program. Terminating.\n");    	
		MPI_Abort(MPI_COMM_WORLD, rc);
    }

	/* get parameters */
	AskParams(&options, argc, argv);              /* ************************* */

	initVariables(&arguments, &results, &options);           /* ******************************************* */

	allocateMatrices(&arguments);        /*  get and initialize variables and matrices  */
	initMatrices(&arguments, &options);            /* ******************************************* */

	gettimeofday(&start_time, NULL);                   /*  start timer         */
	calculate(&arguments, &results, &options);                                      /*  solve the equation  */
	gettimeofday(&comp_time, NULL);                   /*  stop timer          */

	displayStatistics(&arguments, &results, &options);
	DisplayMatrix(&arguments, &results, &options);
	//DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)

	freeMatrices(&arguments);                                       /*  free memory     */
	MPI_Finalize();
	return 0;
}
