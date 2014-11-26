/* HLR WS2014 Übungsblatt 6 Aufgabe 1 HuHufschmidtTimm */

#include <stdio.h>
#include <mpi.h>
#include <sys/time.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int main (int argc, char **argv)
{	
	char str[42];
	int size, rank, rc, dest, source, tag;  
	int i;
	long every_time,min_time;
	int random[100];
	MPI_Status Stat;
	/* mpi starts*/
	rc = MPI_Init (&argc, &argv);	
	if (rc != MPI_SUCCESS) 
	{
		printf ("Error starting MPI program. Terminating.\n");    	
		MPI_Abort(MPI_COMM_WORLD, rc);
    }
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	for (i= 0 ; i < 100 ; i++)
	{
		random[i] = rand();
	}
	/* Rang 0*/
	if (rank == 0)
	{
		every_time = 999999; //make the maximum value.
		for(i=1; i<size; i++)
		{	
			source = i;
			tag = i;
			MPI_Recv(&str, 42, MPI_CHAR, source, tag, MPI_COMM_WORLD, &Stat);
			printf("%s\n", str);			
		}
	}
	/*Rang 1 to n*/
	else
	{	 
		struct timeval time;	
		struct tm *tm;
		char hostname[10];
		gettimeofday(&time, NULL);
		tm = localtime(&time.tv_sec);
		gethostname(hostname, sizeof(hostname));
		every_time = time.tv_usec;
		sprintf(str,"%s: %d-%d-%d %d:%02d:%02d.%ld", hostname, tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday, tm->tm_hour, tm->tm_min, tm->tm_sec, time.tv_usec);		
		dest = 0;
		tag = rank;
		MPI_Send(&str, 42, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
	}	
	
	MPI_Reduce( &every_time, &min_time, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		printf("%ld\n", min_time);
	}
	/*Barrier*/
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Rang %d beendet jetzt!\n", rank);
	
	MPI_Finalize();
	return 0;
}
