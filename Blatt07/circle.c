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
//circle (int* buf, int rank, int size, int* N_per_rank)
circle (int* buf, int rank, int size, int N)
{
	//todo
	int tag = 1234;
	MPI_Status Stat[2];
	MPI_Request Request[2];
	if (rank == 0)
	{
////		printf("Rank 0 Sending");
		//MPI_Isend(buf, N_per_rank[rank], MPI_INT, rank+1, tag, MPI_COMM_WORLD, &Request[0]);
		//MPI_Irecv(buf, N_per_rank[size-1], MPI_INT, size-1, tag, MPI_COMM_WORLD, &Request[1]);
	//}
	//else if (rank == size -1)
	//{
////		printf("last Rank Sending");
		//MPI_Isend(buf, N_per_rank[size-1], MPI_INT, 0, tag, MPI_COMM_WORLD, &Request[0]);
		//MPI_Irecv(buf, N_per_rank[rank-1], MPI_INT, rank-1, tag, MPI_COMM_WORLD, &Request[1]);
	//}
	//else 
	//{	
////		printf("%d, Sending",rank);
		//MPI_Isend(buf, N_per_rank[rank], MPI_INT, rank+1, tag, MPI_COMM_WORLD, &Request[0]);
		//MPI_Irecv(buf, N_per_rank[rank-1], MPI_INT, rank-1, tag, MPI_COMM_WORLD, &Request[1]);
	//}
	//		printf("Rank 0 Sending");
		MPI_Isend(buf, N/size +1 , MPI_INT, rank+1, tag, MPI_COMM_WORLD, &Request[0]);
		MPI_Irecv(buf,  N/size +1 , MPI_INT, size-1, tag, MPI_COMM_WORLD, &Request[1]);
	}
	else if (rank == size -1)
	{
//		printf("last Rank Sending");
		MPI_Isend(buf,  N/size +1 , MPI_INT, 0, tag, MPI_COMM_WORLD, &Request[0]);
		MPI_Irecv(buf,  N/size +1 , MPI_INT, rank-1, tag, MPI_COMM_WORLD, &Request[1]);
	}
	else 
	{	
//		printf("%d, Sending",rank);
		MPI_Isend(buf,  N/size +1, MPI_INT, rank+1, tag, MPI_COMM_WORLD, &Request[0]);
		MPI_Irecv(buf, N/size +1 , MPI_INT, rank-1, tag, MPI_COMM_WORLD, &Request[1]);
	}
//	printf("%d,Waiting",rank);
	//MPI_Waitall(2,&Request,&Stat);
	MPI_Waitall(2,Request,Stat);
	return buf;
}

int
main (int argc, char** argv)
{
	char arg[256];
	int N;	 
	int size, rank, rc, i, rest;
	int* buf;
	int* N_per_rank;	//elements number of each process.
	int do_cycle= 1;
	int a;	//store the first element of process 0.
	int tag=999, tag1=99; //distinguish lines in Vampir from different tag. 
	int flag;
	int displacement=0; //for Outputs after circle. To calculate after circle how many displacements does a array take between processes.
	int origin_rank; // After circle, the present array comes from which rank at initial time.
	
	//int msglen;
	MPI_Status Stat;
	
	if (argc < 2)
	{
		printf("Arguments error\n");
		return EXIT_FAILURE;
	}
	sscanf(argv[1], "%s", arg);
	//array length
	N = atoi(arg);
	
	/* mpi starts*/
	rc = MPI_Init (&argc, &argv);	
	if (rc != MPI_SUCCESS) 
	{
		printf ("Error starting MPI program. Terminating.\n");    	
		MPI_Abort(MPI_COMM_WORLD, rc);
    	}
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	
  	/*array division scheme*/
	
	if (N < size) 
	{
	printf("Array length shorter than number of processes. This is forbidden. Try using the program with a greater array length \n");
	return -2;
	}
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
  	//printf("%d,Cycle here",rank);
	while(do_cycle)
	{
//		printf("%d, calling circle",rank);
		circle(buf,rank,size,N);
		//MPI_Barrier(MPI_COMM_WORLD);
		if(rank==size-1)
		{
			if (buf[0]==a)
			{			
				do_cycle = 0;	
			//}
				//printf("Sending Abort");
				//MPI_Request *Request2 = malloc(sizeof(MPI_Request *)*size);
				MPI_Request *Request2;
				MPI_Status *Stat2;
				Request2 = (MPI_Request *) malloc (sizeof(MPI_Request) * (size-1));
				//MPI_Status *Stat2 = malloc(sizeof(MPI_Status *) * size);
				Stat2=(MPI_Status *) malloc(sizeof(MPI_Status) * (size-1));
				
				//predefine the Requests to prevent segfaults
				for (i=0; i<size-1; ++i)
				{
				 Request2[i]=MPI_REQUEST_NULL;
				} 
				for (i=0; i<size-1; i++)
			//{
			//	MPI_Send(&do_cycle, 1, MPI_INT, i, tag1, MPI_COMM_WORLD);
			//}
				{
					MPI_Isend(&do_cycle, 1, MPI_INT, i, tag1, MPI_COMM_WORLD, &Request2[i]);
					//MPI_Wait(&Request2[i],&Stat2[i]);
				}
				MPI_Waitall(size-1,Request2,Stat2);
			}	
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(rank!=size-1)
		//else
		//{
		//	//MPI_Iprobe(size-1,tag1,MPI_COMM_WORLD,&flag,&Stat);
		//	//if (flag) {
		//	  	//MPI_Get_count(&Stat, MPI_INT, &msglen);
		//		MPI_Recv(&do_cycle, 1, MPI_INT, size-1, tag1, MPI_COMM_WORLD, &Stat);
		//	//}
		//}
		{
			MPI_Status Stat2;
			MPI_Request Request2;
//			printf("%d,Test Abrot Msg",rank);
			MPI_Iprobe(size-1,tag1,MPI_COMM_WORLD,&flag,&Stat2);
			if (flag)
				{
			  	//MPI_Get_count(&Stat, MPI_INT, &msglen);
				//printf("%d, Recive MSG",rank);
				MPI_Irecv(&do_cycle, 1, MPI_INT, size-1, tag1, MPI_COMM_WORLD, &Request2);
				MPI_Wait(&Request2,&Stat2);
				}
		}
		displacement++; //for Outputs after circle.
	
		//MPI_Barrier(MPI_COMM_WORLD);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	printf("\nAFTER\n");

	origin_rank = rank-displacement;
	if(origin_rank < 0)
	{
		origin_rank = origin_rank +size;
	}
	//printf("i_dis:%dfrom %d\n", i_dis,rank);
	//printf("N_per_rank[rank]%d\n", N_per_rank[rank]);
	//printf("N_per_rank[rank+1]%d\n", N_per_rank[rank+1]);
	//printf("N_per_rank[i_dis]%d\n", N_per_rank[i_dis]);
	if (origin_rank < rest)
	{
		for (int j = 0; j < N/size + 1; j++)
		{
			printf ("rank %d: %d\n", rank, buf[j]);
		}
	}
	else
	{
		for (int j = 0; j < N/size; j++)
		{
			printf ("rank %d: %d\n", rank, buf[j]);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	free(N_per_rank);
	MPI_Finalize();
	return EXIT_SUCCESS;
}
