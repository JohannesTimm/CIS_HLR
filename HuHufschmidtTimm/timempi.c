#include mpi.h
#include <stdio.h>
#include <stdlib.h>

#define STRINGLENTH 100
int main(int argc, char** argv)) {
	//Init MPI
	MPI_Init(&argc,&argv);
	// get size
	int size_of_world;
	MPI_Comm_size(MPI_COMM_WORLD, &size_of_world);
	//get rank
	int rank_in_world;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_in_world);
	
	char string_buffer[STRINGLENTH];
	char* read_buffer =NULL;
	if (rank_in_world == 0 )
		{
		read_buffer=malloc(size_of_world * STRINGLENGTH *sizeof(char));
		}
	char name_of_host[STRINGLENTH];
	int name_length=STRINGLENTH;
	//hostname
	if (gethostname(name_of_host,name_length) !=0)
		{
		printf("Error while gathering the hostname");
		exit(1);
		}
	//time
	struct timeval time;
	if (gettimeofday(&time,NULL) !=0) 
		{
		printf("Error while gahtering the time");
		exit(2);
		}
	//Output
	time_t current_time=time.tv_sec;
	char time_buffer[30] //size choosen by rolling dice
	strftime(time_buffer,30, "%Y-%m-%d %T.",localtime(&current_time));
	sprintf(string_buffer, "%s: %s%li", name_of_host, time_buffer, time.tv_usec);
	//Gather that stuff
	int rc;
	rc=MPI_Gather(string_buffer, LEN,MPI_CHAR, read_buffer, LEN,MPI_CHAR,0, MPI_COMM_WORLD);
	if (rc!= MPI_SUCCESS)
		{
		printf("Error while MPI_Gather, Error Code is : %d", rc);
		exit(3);
		}
	//print outputs
	if (rank_in_world == 0) 
		{
		for (int i=0;i< size_of_world; ++i)
			{
			printf(".*s\n", LEN,read_buffer + LEN *i);
			}
		}
	//Calculate Microseconds
	int microsec;
	microsec= time.tv_usec;
	int * read_buffer_microsec;
	if (rank_in_world == 0) 
		{
		read_buffer_microsec=malloc(1 * sizeof(int);
		}
	//MPI_REDUCE
	int rc;
	rc=MPI_Reduce(&microsec,read_buffer_microsec,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
	if (rc!= MPI_SUCCESS)
		{
		printf("Error while MPI_Reduce, Error Code is : %d", rc);
		exit(4);
		}
	//print them
	if (rank_in_world == 0) 
		{
		printf("d\n", read_buffer_microsec);
		}
	
	printf("Process from host %d, Rang %b beendt jetzt!\n", name_of_host,rank_in_world);
	
	MPI_Finalize();