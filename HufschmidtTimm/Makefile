CC=gcc
CFLAGS=-std=gnu99 -Wall -Wextra -pedantic
LD=$(CC)
LDFLAGS=-lm

all: map pointer

pointer: pointer.o

map: map.o

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f map pointer
	rm -f *.o *~
