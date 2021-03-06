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
#include <omp.h>
#include "partdiff-par-hybrid.h"

struct calculation_arguments
{
	//uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	//uint64_t  num_matrices;   /* number of matrices                             */
	int  N;              /* number of spaces between lines (lines=N+1)     */
	int  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
	int N_local;
	int rank;
	int size;
	int from;
	int to;
	//uint64_t N_local;
	//uint64_t rank;
	//uint64_t size;
	//uint64_t from;
	//uint64_t to;
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision_local; /* actual precision of all slaves in iteration    */
	double    stat_precision;
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
	//uint64_t rest;
	//uint64_t const N = arguments->N;
	
	MPI_Comm_rank (MPI_COMM_WORLD, &arguments -> rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &arguments -> size);
	
	int const size = arguments -> size;
	int const rank = arguments -> rank;	
	//uint64_t const size = arguments -> size;
	//uint64_t const rank = arguments -> rank;	
	
	arguments->N = (options->interlines * 8) + 9 - 1;
	//uint64_t const N = arguments->N;
	int const N = arguments->N;
	//N = (options->interlines * 8) + 9 - 1;
	rest = (N + 1 -2) % size;
	
	//printf("Global Matrix Size is %d the rest is %d \n",(int)arguments->N + 1,rest); 
	//printf("Global Matrix Size is %d the rest is %d \n",N + 1,rest); 	
	if(rank < rest)
	{
		arguments -> N_local = (N + 1 -2) / size + 1;
		arguments -> from = rank * ((N - 1)/size + 1) + 1;
			//arguments -> from = rank * ((N - 1)/size + 1)+2;
		arguments -> to = arguments -> from + (N-1)/size;
	}
	else 
	{
		arguments -> N_local = (N + 1 -2) / size;
		arguments -> from = rank * ((N - 1)/size) + rest + 1;
		arguments -> to = arguments -> from + (N-1)/size - 1;
	}
	//printf("Process %d has local Matrix size of %d and is responsible for Global lines %d to %d \n",rank, arguments->N_local,arguments->from,arguments->to); 
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
	//uint64_t i;
	int i;
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
allocateMemory (size_t sizet)
{
	void *p;

	if ((p = malloc(sizet)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes)\n", sizet);
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
	//uint64_t i, j;

	//uint64_t const N = arguments->N;
	//uint64_t const N_local = arguments -> N_local;
	int i, j;

	int const N = arguments->N;
	int const N_local = arguments -> N_local;

	arguments->M = allocateMemory(arguments->num_matrices * (N_local + 2) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; ++i)
	{
		arguments->Matrix[i] = allocateMemory((N_local + 2) * sizeof(double*));

		for (j = 0; j < N_local + 2 ; ++j)
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
	//uint64_t g, i, j;                                /*  local variables for loops   */

	//uint64_t const N = arguments -> N;
	//uint64_t const N_local =arguments -> N_local;
	int g, i, j;                                /*  local variables for loops   */

	int const N = arguments -> N;
	int const N_local =arguments -> N_local;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;
	int from = arguments -> from;
	int const rank = arguments -> rank;
	int const size = arguments -> size;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; ++g)
	{
		//for (i = 0; i < N_local + 2; ++i)
		for (i = 0; i <= N_local + 1; ++i)
		{
			for (j = 0; j <= N; ++j)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; ++g)
		{
			//for (i = 0; i < N_local + 2; i++)
			for (i = 0; i < N_local + 1; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i) - (from - 1) * h;
				Matrix[g][i][N] = h * i + (from - 1) * h;
			}
			
			if (rank == 0)	
			{	
				for (i = 0; i <= N; ++i)
				{
					Matrix[g][0][i] = 1.0 - (h * i);
				}
			}
			if (rank + 1 == size)
			{
				for (i = 0; i <= N; ++i)
				{
					Matrix[g][N_local + 1][i] = h * i;
					//Matrix[g][N_local][i] = h * i;
				}
			}
			//for (i = 0; i <= N; ++i)
			//{
				//Matrix[g][N_local + 1][i] = h * i;
				////Matrix[g][N_local][i] = h * i;
				//Matrix[g][0][N] = 0.0;
			//}
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* For single process single thread only - compatibility code		    */
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

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/*     */
/* ************************************************************************ */
static
void
calculate_MPI_Jacobi (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options)
{
	int i, j;                                   /* local variables for loops  */
	int k;
	int m1, m2;                                 /* used as indices for old and new matrices       */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration                  */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */
	double globalmaxresiduum;
	int const N_local =arguments->N_local;
	int const N = arguments->N;
	MPI_Status status;
	int const Max_Threads=options->number;
 	int from = arguments -> from;
 	int to = arguments -> to;
 	int const size = arguments -> size;
  	int const rank = arguments -> rank;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;
	int error_code;
	int tag_send, tag_recv;
	int thread_level, thread_is_main;
	omp_set_dynamic(0); //disable dynamic teams to enforce thread number
	omp_set_num_threads(Max_Threads);
	
	MPI_Query_thread(&thread_level);
	//MPI_Is_thread_main(&thread_is_main);

	
	//printf("Rank %d, N_local %d, Interval %d, From %d, To %d \n",rank,N_local,(to-from),from, to);

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
		int rest;
		int i_start, i_end;
		int num_threads;
		int my_thread;
		int width;
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];
// 		double* fpisin_i;
// 		fpisin_i = (double*)malloc(N * sizeof(double));
		maxresiduum = 0;
//		#pragma omp parallel shared(Matrix_In,Matrix_Out) private(my_thread,i_start,i_end,thread_is_main) reduction(max:maxresiduum)
//		#pragma omp parallel shared(Matrix_In,Matrix_Out), private(j,star,residuum), reduction(max:maxresiduum)
		#pragma omp parallel private(j,star,residuum, i_start,i_end,thread_is_main)
		{
			
			
			num_threads = omp_get_num_threads();
			my_thread = omp_get_thread_num();
// 			width = (int) (N_local) / num_threads;
// 			rest = (N_local) % num_threads;
			MPI_Is_thread_main(&thread_is_main);
			//if(my_thread == num_threads-1)
		
			//{
			//	i_start = 1 + my_thread * width;
			//	i_end = N_local;
			//}
			//else
			//{
			//	i_start = 1 + my_thread * width;
			//	i_end = i_start + width;
			//}
// 			if(my_thread < rest)
// 				{
// 				i_start=my_thread*((N_local)/num_threads +1)+1;
// 				i_end=i_start+(N_local)/size;
// 				}
// 			else
// 				{
// 				i_start=my_thread*((N_local)/num_threads) +rest+1;
// 				i_end=i_start+(N_local)/size-1;
// 				}	
// 			if (options->inf_func == FUNC_FPISIN)
// 			{
// 				#pragma omp for firstprivate(fpisin,pih) //,i_start,i_end)
// 				for (i = i_start; i < i_end; i++)			
// 				{ 	
// 					fpisin_i[i] = 0.0;
// 					fpisin_i[i]= fpisin * sin(pih * ((double)i + from - 1));
// 				}
// 			}
			//printf("Here is Thread %d of %d from Rank %d Master %d istart %d iend %d \n",my_thread,num_threads,rank,thread_is_main,i_start,i_end);
			/* over all rows */
			double fpisin_i = 0.0;
// 			#pragma omp for private(j,star,residuum) firstprivate(pih) reduction(max:maxresiduum) collapse(2)
			#pragma omp for private(j,star,residuum,fpisin_i) firstprivate(pih,fpisin) reduction(max:maxresiduum)
			for (i = 1; i <= N_local; ++i)
			//for (i = i_start; i <= i_end; ++i)
			{	
				//printf("Thread %d of %d from Process %d Calculating Row %d \n",my_thread,num_threads,rank,i);
				if (options->inf_func == FUNC_FPISIN)
				{
					fpisin_i = fpisin * sin(pih * ((double)i + from - 1));
				}

				/* over all columns */
				for (j = 1; j < N; ++j)
				{
					star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);	
	
					if (options->inf_func == FUNC_FPISIN)
					{
						star += fpisin_i * sin(pih * (double)j);
					}
// 					if (options->inf_func == FUNC_FPISIN)
// 					{
// 						star += fpisin_i[i] * sin(pih * (double)j);
// 					}
	
					if (options->termination == TERM_PREC || term_iteration == 1)
					{
						residuum = Matrix_In[i][j] - star;
						residuum = (residuum < 0) ? -residuum : residuum;
						maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
					}
	
					Matrix_Out[i][j] = star;
				}
			}
			#pragma omp barrier
			#pragma omp master
			{	
				MPI_Query_thread(&thread_level);
				MPI_Is_thread_main(&thread_is_main);
				if ((thread_level > MPI_THREAD_FUNNELED ) || (thread_level == MPI_THREAD_FUNNELED && thread_is_main))
				{
					//Communicate with the Other Processes to exchange Halo lines
					if (rank>0)
					{
						//tag_send=rank+10;
						//tag_recv=rank-1+20;
						tag_send=rank;
						tag_recv=rank-1;
						error_code=MPI_Sendrecv(Matrix_Out[1],N + 1,MPI_DOUBLE, rank -1, tag_send,
									Matrix_Out[0],N + 1,MPI_DOUBLE, rank -1, tag_recv,
									MPI_COMM_WORLD,&status);
						if (error_code!=MPI_SUCCESS)
						{
							printf("Error in SendRecv");
						}
					}
					if (rank<size -1)
						{
						tag_send=rank;
						tag_recv=rank+1;
						/*error_code=MPI_Sendrecv(Matrix_Out[N_local + 2 - 2],N + 1,MPI_DOUBLE, rank +1, tag_send,
									Matrix_Out[N_local + 2 - 1],N + 1,MPI_DOUBLE, rank +1, tag_recv,
										MPI_COMM_WORLD,&status);
						*/			
						error_code=MPI_Sendrecv(Matrix_Out[N_local],N + 1,MPI_DOUBLE, rank +1, tag_send,
								Matrix_Out[N_local+1],N + 1,MPI_DOUBLE, rank +1, tag_recv,
									MPI_COMM_WORLD,&status);			
						if (error_code!=MPI_SUCCESS)
						{
								printf("Error in SendRecv");
						}
					}
					//Find the global maximum of the residuum (later decide if this is low enough)
					//Allreduce has no error handling
					MPI_Allreduce(&maxresiduum,&globalmaxresiduum,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			
		
					/* exchange m1 and m2 */
					k = m1;
					m1 = m2;
					m2 = k;
			
					results->stat_iteration++;
					results->stat_precision_local = maxresiduum;
					results->stat_precision = globalmaxresiduum;
				
					MPI_Barrier(MPI_COMM_WORLD);
				}
				else
				{
					printf("Error: Master Thread of Process %d cannot communicate. Why?",rank);
				}
				
				
				/* check for stopping calculation, depending on termination method */
				if (options->termination == TERM_PREC)
				{
					if (globalmaxresiduum < options->term_precision)
					{
					term_iteration = 0;
						
					}
				}
				else if (options->termination == TERM_ITER)
				{
					
					term_iteration--;
				}
			}
			#pragma omp barrier	
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
  double* result_vec;
  MPI_Status status;
  
  int from = arguments -> from;
  int to = arguments -> to;
  
  int const size = arguments -> size;
  int const rank = arguments -> rank;
  result_vec=malloc(elements*sizeof(double));
  /* first line belongs to rank 0 */
  if (rank == 0)
    from--;

  /* last line belongs to rank size - 1 */
  if (rank + 1 == size)
    to++;
	
  //printf("Rank %d, responsible for lines from %d to %d \n", rank, from, to);	
  
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
        //printf("Recv Line, %d , with tag , %d",line,42+y); 
        //MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
        MPI_Recv(result_vec, elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
        //printf("Got Line");
      }
      else
      {
      //printf("Rank %d:",rank);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        /* if the line belongs to this process, send it to rank 0
         * (line - from + 1) is used to calculate the correct local address */
        //printf("Send Line %d, from %d, with tag %d", line,rank,42+y);
        //printf("Rank %d Send Line %d Local %d from %d to %d \n", rank, line, (line-from+1),from,to);
        MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
        //printf("Line Send");
      }
    }

    if (rank == 0)
    {
      //printf("L %d:",line);
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
          //printf("%7.4f", Matrix[0][col]);
          printf("%7.4f", result_vec[col]);
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
	int thread_level_required=MPI_THREAD_MULTIPLE;

	int thread_level;
	/* mpi starts*/
	rc= MPI_Init_thread(&argc,&argv,thread_level_required,&thread_level);

	if (rc != MPI_SUCCESS) 
	{
		printf ("Error starting MPI program. Terminating.\n");    	
		MPI_Abort(MPI_COMM_WORLD, rc);
    	}
    	MPI_Query_thread(&thread_level);
	if ((thread_level > MPI_THREAD_FUNNELED ) || (thread_level == MPI_THREAD_FUNNELED ))
		{
			//printf("Using Threads\n");
		}
	else
		{
			printf("Error : Thread Level provided by Libary is too low. Programm will not work\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
    	}
	
	/* get parameters */
	AskParams(&options, argc, argv);              /* ************************* */

	initVariables(&arguments, &results, &options);           /* ******************************************* */

	allocateMatrices(&arguments);        /*  get and initialize variables and matrices  */
	initMatrices(&arguments, &options);            /* ******************************************* */
	MPI_Barrier(MPI_COMM_WORLD); //all should start at nearly the same time
	gettimeofday(&start_time, NULL);                   /*  start timer         */
	if (options.method == METH_JACOBI)
	{
		//calculate(&arguments, &results, &options);                                     /*  solve the equation  */
		calculate_MPI_Jacobi(&arguments, &results, &options);
	}
	else
	{	
		printf("The Gauss-Seidel Method is not yet implementet as a parallel program. Computing ressources are Wasted");
		calculate(&arguments, &results, &options);     
	}
	gettimeofday(&comp_time, NULL);                   /*  stop timer          */
	MPI_Barrier(MPI_COMM_WORLD);
	if (arguments.rank==0) //just have the statistics once and not for every process
	{
		displayStatistics(&arguments, &results, &options);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	DisplayMatrix(&arguments, &results, &options);
	//DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
	MPI_Barrier(MPI_COMM_WORLD);
	freeMatrices(&arguments);                                       /*  free memory     */
	MPI_Finalize();
	return 0;
}
