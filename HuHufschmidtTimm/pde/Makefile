# Common definitions
CC = gcc

# Compiler flags, paths and libraries
# CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

PROGS = partdiff-openmp partdiff-openmp-zeilen partdiff-openmp-spalten partdiff-openmp-element

OBJS = askparams.o

# Targets ...
all: $(PROGS)

partdiff-openmp: $(OBJS) partdiff-openmp.o 
	$(CC) $(LFLAGS) -o $@ $(OBJS) partdiff-openmp.o  $(LIBS)

partdiff-openmp.o: partdiff-openmp.c
	$(CC) -c $(CFLAGS) partdiff-openmp.c

partdiff-openmp-zeilen: $(OBJS) partdiff-openmp-zeilen.o 
	$(CC) $(LFLAGS) -o $@ $(OBJS) partdiff-openmp-zeilen.o  $(LIBS)

partdiff-openmp-zeilen.o: partdiff-openmp.c
	$(CC) -c $(CFLAGS) -DAUFTEILUNG=1 -o partdiff-openmp-zeilen.o  partdiff-openmp.c

partdiff-openmp-spalten: $(OBJS) partdiff-openmp-spalten.o 
	$(CC) $(LFLAGS) -o $@ $(OBJS) partdiff-openmp-spalten.o  $(LIBS)

partdiff-openmp-spalten.o: partdiff-openmp.c
	$(CC) -c $(CFLAGS) -DAUFTEILUNG=2 -o partdiff-openmp-spalten.o  partdiff-openmp.c

partdiff-openmp-element: $(OBJS) partdiff-openmp-element.o 
	$(CC) $(LFLAGS) -o $@ $(OBJS) partdiff-openmp-element.o  $(LIBS)

partdiff-openmp-element.o: partdiff-openmp.c
	$(CC) -c $(CFLAGS) -DAUFTEILUNG=3 -o partdiff-openmp-element.o  partdiff-openmp.c

askparams.o: askparams.c
	$(CC) -c $(CFLAGS) askparams.c

clean:
	$(RM) *.o *~ $(PROGS)

