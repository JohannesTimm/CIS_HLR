#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int*
init (int N)
{
  //todo
  int* buf = malloc(sizeof(int) * N);

  srand(time(NULL));

  for (int i = 0; i < N; i++)
  {
    buf[i] = rand() % 25; //do not modify %25
  }

  return buf;
}

int*
circle (int* pbuf)
{
  //todo
  
  return buf;
}

int
main (int argc, char** argv)
{
  //char arg[256];
  int rank;
  int N;
  int* buf;
  int rc,size;

  //if (argc < 2)
  //{
    //printf("Arguments error\n");
    //return EXIT_FAILURE;
  //}

  //sscanf(argv[1], "%s", arg);

  //array length
  //N = atoi(arg);
	N= 3;
		rc = MPI_Init (&argc, &argv);	
	if (rc != MPI_SUCCESS) 
	{
		printf ("Error starting MPI program. Terminating.\n");    	
		MPI_Abort(MPI_COMM_WORLD, rc);
    }
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	buf = init(N);
  
  //************************************
	int* displs;
	int* scounts;
	int* pbuf;
	int i, rest;
    displs = (int *)malloc(size*sizeof(int)); 
    scounts = (int *)malloc(size*sizeof(int));
	//sum = (int *)malloc(size*sizeof(int));
	if ((N % size) == 0)
	{	
		for (i=0; i<size; i++)
		{
			scounts[i] = (int)(N / size);
			//if (rank ==0)
			//{
			//printf("%d\n", scounts[i]);
			//}
		}
	}
 	else 
	{	
		rest = N % size;
		for (i=0; i<rest; i++)	
		{
			scounts[i] = (int)(N / size) + 1;		
		}
		for (i=rest; i<size; i++)
		{	
			scounts[i] = (int)(N / size) ;
		}
	}
	pbuf = (int *)malloc(sizeof(int)*(N/size + 1));
    //sum[0] = 0;
    displs[0] = 0;
	for (i=1; i<size; i++)
	{ 
        displs[i] = displs[i-1] + scounts[i-1];    
        //if (rank == 0)
        //{
			//for (i=0; i<size; i++)
			//{
				//printf("%d\n", displs[i]);
			//}
		//}
        //sum[i] = sum[i] + scounts[i]  
	} 
	
    /* Create datatype for the column we are receiving 
     */ 
    //MPI_Type_vector( 100-myrank, 1, 150, MPI_INT, &rtype); 
    //MPI_Type_commit( &rtype ); 
    //rptr = &recvarray[0][myrank]; 
	MPI_Scatterv( buf, scounts, displs, MPI_INT, pbuf, N/size + 1, MPI_INT, 0, MPI_COMM_WORLD); 

//*********************************************************

  //todo myrank
  //rank = 0;

  printf("\nBEFORE\n");

  for (int i = 0; i < scounts[rank]; i++)
  {
    printf ("rank %d: %d\n", rank, pbuf[i]);
  }

  circle(pbuf);

  printf("\nAFTER\n");

  for (int j = 0; j < scounts[rank]; j++)
  {
    printf ("rank %d: %d\n", rank, pbuf[j]);
  }
		MPI_Finalize();
	 return EXIT_SUCCESS;
}
