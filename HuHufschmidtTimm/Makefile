# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJSMPI = partdiff-mpi.o askparams.o

OBJHYBRID = partdiff-hybrid.o askparams.o
# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-mpi partdiff-hybrid

partdiff-mpi: $(OBJSMPI) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJSMPI) $(LIBS)
partdiff-hybrid: $(OBJHYBRID) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJHYBRID) $(LIBS)
clean:
	$(RM) *.o *~

partdiff-mpi.o: partdiff-mpi.c Makefile
partdiff-hybrid.o: partdiff-hybrid.c Makefile
askparams.o: askparams.c Makefile
