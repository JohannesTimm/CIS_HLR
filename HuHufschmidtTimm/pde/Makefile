# Common definitions
CC = gcc

# Compiler flags, paths and libraries
# CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS = partdiff-openmp.o askparams.o

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-openmp

partdiff-openmp: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	$(RM) *.o *~

<<<<<<< HEAD
partdiff-openmp.o: partdiff-openmp.c Makefile
=======
partdiff-seq.o: partdiff-openmp.c Makefile
>>>>>>> 2bb028815b10846b98b843253d55f53967baa410

askparams.o: askparams.c Makefile
