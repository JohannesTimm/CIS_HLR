# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJSMPI = partdiff-par.o askparams.o

OBJHYBRID = partdiff-par-hybrid.o askparams.o
# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-par
#partdiff-par-hybrid

partdiff-par: $(OBJSMPI) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJSMPI) $(LIBS)
partdiff-par-hybrid: $(OBJHYBRID) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJHYBRID) $(LIBS)
clean:
	$(RM) *.o *~

partdiff-par.o: partdiff-par.c Makefile
partdiff-par-hybrid.o: partdiff-par-hybrid.c Makefile
askparams.o: askparams.c Makefile
