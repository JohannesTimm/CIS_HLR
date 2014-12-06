VAMPCC="vtcc -vt:cc mpicc -std=c99 -ggdb -Wall -Wextra"
NUMPROC="3"
circle: circle.c Makefile
	mpicc -std=c99 -g -ggdb -Wall -Wextra -fbounds-check circle.c -o circle
vampir: circle.c
	#$(VAMPCC) circle.c -o circle
	vtcc -vt:cc mpicc -std=c99 -ggdb -Wall -Wextra circle.c -o circle
	mpiexec -np $(NUMPROC) ./circle $(NUMPROC)
	vampir circle.otf
clean:
	rm -f circle
