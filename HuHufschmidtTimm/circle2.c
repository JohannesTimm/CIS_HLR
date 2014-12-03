#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int*
init (int N, int rank, int size, int rest, int* N_per_rank)
{
	//todo
	int* buf;
	if ((N % size)==0)
	{			
		buf = malloc(sizeof(int) * N_per_rank[rank]);
		srand(time(NULL));

		for (int i = 0; i < N_per_rank[rank]; i++)
		{
			buf[i] = rand() % 25; //do not modify %25
		}
	}
	else 
	{	
		if (rank < rest)
		{
			buf = malloc(sizeof(int) * N_per_rank[rank]);
			srand(time(NULL));
			for (int i = 0; i < N_per_rank[rank]; i++)
			{
				buf[i] = rand() % 25; //do not modify %25
			}
		}
		else
		{	
			buf = malloc(sizeof(int) * N_per_rank[rank]);
			srand(time(NULL));
			for (int i = 0; i < N_per_rank[rank] ; i++)
			{
				buf[i] = rand() % 25; //do not modify %25
			}
		}
	}
  return buf;
}

int*
//circle (int* buf, int rank, int size, int* N_per_rank)
circle(int* buf)
{
  //todo
	//int tag = rank + 10;
	//MPI_Status Stat;
	//MPI_Request Request;
	//if (rank == 0)
	//{
		//MPI_Isend(buf, N_per_rank[rank], MPI_INT, rank+1, tag, MPI_COMM_WORLD, &Request);
		//MPI_Irecv(buf, N_per_rank[size-1], MPI_INT, size-1, tag, MPI_COMM_WORLD, &Request);
	//}
	//else if (rank == size -1)
	//{
		//MPI_Isend(buf, N_per_rank[0], MPI_INT, 0, tag, MPI_COMM_WORLD, &Request);
		//MPI_Irecv(buf, N_per_rank[size-1], MPI_INT, size-1, tag, MPI_COMM_WORLD, &Request);
	//}
	////else 
	////{
		////MPI_Isend(buf, N_per_rank[rank], MPI_INT, rank+1, tag, MPI_COMM_WORLD, &Request);
		////MPI_Irecv(buf, N_per_rank[size-1], MPI_INT, size-1, tag, MPI_COMM_WORLD, &Request);
	////}
	////MPI_Wait(&Request, &Stat);
	return buf;
}

int
main (int argc, char** argv)
{
	char arg[256];
	int N;
	int size;
	int rank;
	int rc;
	int* buf;
	int* N_per_rank;
	int rest;

	//if (argc < 2)
	//{
		//printf("Arguments error\n");
		//return EXIT_FAILURE;
	//}

	//sscanf(argv[1], "%s", arg);

	//array length
	//N = atoi(arg);
	N = 5;
	
	/* mpi starts*/
	rc = MPI_Init (&argc, &argv);	
	if (rc != MPI_SUCCESS) 
	{
		printf ("Error starting MPI program. Terminating.\n");    	
		MPI_Abort(MPI_COMM_WORLD, rc);
    }
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);

  //todo myrank
  //rank = 0;
	N_per_rank = malloc(sizeof(int) * size);
	if ((N % size) == 0)
	{
		N_per_rank[rank] = (int)(N / size);
	}
 	else 
	{	
		rest = N % size;
		if (rank < rest)
		{
			N_per_rank[rank] = (int)(N / size) + 1;		
		}
		else
		{	
			N_per_rank[rank] = (int)(N / size) ;
		}
	}
	buf = init(N,rank,size,rest,N_per_rank);
	printf("\nBEFORE\n");

  for (int i = 0; i < N_per_rank[rank]; i++)
  {
    printf ("rank %d: %d\n", rank, buf[i]);
  }

  //circle(buf,rank,size,N_per_rank);

  printf("\nAFTER\n");

  for (int j = 0; j < N_per_rank[rank]; j++)
  {
    printf ("rank %d: %d\n", rank, buf[j]);
  }
	MPI_Finalize();
  return EXIT_SUCCESS;
}
