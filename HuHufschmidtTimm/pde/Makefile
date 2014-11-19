# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb
EXTRACC = -pthread
LFLAGS = $(CFLAGS)
LIBS   = -lm
EXTRALD = -pthread
PROGS = partdiff-posix partdiff-seq
OBJS = askparams.o


# Targets ...
all: $(PROGS)

partdiff-seq: $(OBJS) partdiff-seq.o
	$(CC) $(CFLAGS) -o partdiff-seq  $(OBJS) partdiff-seq.o  $(LIBS)

partdiff-seq.o: partdiff-seq.c
	$(CC) -c $(CFLAGS)  -DMODE=0 -o partdiff-seq.o partdiff-seq.c

partdiff-posix: $(OBJS) partdiff-posix.o
	$(CC) $(LFLAGS) $(EXTRALD)-o $@ $(OBJS) partdiff-posix.o  $(LIBS)

partdiff-posix.o: partdiff-seq.c
	$(CC) -c $(CFLAGS)  $(EXTRACC) -DMODE=1 -o partdiff-posix.o  partdiff-seq.c

askparams.o: askparams.c Makefile
	$(CC) -c $(CFLAGS) askparams.c

clean:
	$(RM) *.o *~ $(PROGS)

