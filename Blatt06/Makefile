CC = mpicc
#CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb
CFLAGS = -std=gnu99 -pedantic -Wall -Wextra -O2 -ggdb
PROGS = timempi timempi2

all: $(PROGS)
timempi: timempi.c
	$(CC) $(CFLAGS) -o timempi timempi.c
timempi2: timempi2.c
	$(CC) $(CFLAGS) -o timempi2 timempi2.c
timempi2-debug:timempi2.c
	$(CC) $(CFLAGS) -DMODE=1 -o timempi2-debug timempi2.c
clean:
	rm $(PROGS) timempi2-debug
