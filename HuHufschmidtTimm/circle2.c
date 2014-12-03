#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int*
init (int N, int rank, int size, int* N_per_rank)
{
	//todo
	int* buf;	
		buf = malloc(sizeof(int) * (N/size)+1);
		srand(time(NULL)+rank);

		for (int i = 0; i < N_per_rank[rank]; i++)
		{
			buf[i] = rand() % 25; //do not modify %25
		}
	return buf;
}

int*
circle (int* buf, int rank, int size, int* N_per_rank)
{
	//todo
	int tag = 1234;
	MPI_Status Stat[2];
	MPI_Request Request[2];
	if (rank == 0)
	{
		MPI_Isend(buf, N_per_rank[rank], MPI_INT, rank+1, tag, MPI_COMM_WORLD, &Request[0]);
		MPI_Irecv(buf, N_per_rank[size-1], MPI_INT, size-1, tag, MPI_COMM_WORLD, &Request[1]);
	}
	else if (rank == size -1)
	{
		MPI_Isend(buf, N_per_rank[size-1], MPI_INT, 0, tag, MPI_COMM_WORLD, &Request[0]);
		MPI_Irecv(buf, N_per_rank[rank-1], MPI_INT, rank-1, tag, MPI_COMM_WORLD, &Request[1]);
	}
	else 
	{	
		MPI_Isend(buf, N_per_rank[rank], MPI_INT, rank+1, tag, MPI_COMM_WORLD, &Request[0]);
		MPI_Irecv(buf, N_per_rank[rank-1], MPI_INT, rank-1, tag, MPI_COMM_WORLD, &Request[1]);
	}
	MPI_Waitall(2,Request,Stat);
	
	return buf;
}

int
main (int argc, char** argv)
{
	//char arg[256];
	int N=9;	 //put the dimension of array.	
	int size, rank, rc, i, rest;
	int* buf;
	int* N_per_rank;	//elements number of each process.
	int do_cycle= 1;
	int a;	//store the first element of process 0.
	int tag=999, tag1=99; //distinguish lines in Vampir from different tag. 
	MPI_Status Stat;
	//if (argc < 2)
	//{
		//printf("Arguments error\n");
		//return EXIT_FAILURE;
	//}
	//sscanf(argv[1], "%s", arg);
	//array length
	//N = atoi(arg);
	
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
  	/*array division scheme*/
	N_per_rank = malloc(sizeof(int) * size);
	if ((N % size) == 0)
	{
		for (i=0 ; i<size; i++)
		{
			N_per_rank[i] = (int)(N / size);
		}
	}
 	else 
	{	
		rest = N % size;
		if (rank < rest)
		{
			for (i=0 ; i<rest; i++) 
			{
				N_per_rank[i] = (int)(N / size) + 1;
			}		
		}
		else
		{	
			for (i=rest ; i<size; i++) 
			{
				N_per_rank[i] = (int)(N / size) ;
			}
		}
	}
	buf = init(N,rank,size,N_per_rank);
	
	/*store the value from first element of process 0 to a*/
	if(rank == 0)
	{
		MPI_Send(&buf[0], 1, MPI_INT, size-1, tag, MPI_COMM_WORLD);
	}
	if(rank == size-1)	
	{
		MPI_Recv(&a, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &Stat);
	}
	
	printf("\nBEFORE\n");

	for (int i = 0; i < N_per_rank[rank]; i++)
	{
		printf ("rank %d: %d\n", rank, buf[i]);
	}
  	/*cycling, when condition is satified, send message from the last process to others processes*/ 
	while(do_cycle)
	{
		circle(buf,rank,size,N_per_rank);
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==size-1)
		{
			if (buf[0]==a)
			{			
				do_cycle = 0;	
			}
			for (i=0; i<size-1; i++)
			{
				MPI_Send(&do_cycle, 1, MPI_INT, i, tag1, MPI_COMM_WORLD);
			}
		}
		else
		{
			MPI_Recv(&do_cycle, 1, MPI_INT, size-1, tag1, MPI_COMM_WORLD, &Stat);
		}		
		MPI_Barrier(MPI_COMM_WORLD);
	}

	printf("\nAFTER\n");

	for (int j = 0; j < N_per_rank[rank]; j++)
	{
		printf ("rank %d: %d\n", rank, buf[j]);
	}
	MPI_Finalize();
	return EXIT_SUCCESS;
}
