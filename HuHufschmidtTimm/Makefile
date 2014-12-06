$VAMPCC="vtcc -vt:cc mpicc -std=c99 -ggdb -Wall -Wextra"
$NUMPROC=3
circle: circle2.c Makefile
	mpicc -std=c99 -g -ggdb -Wall -Wextra -fbounds-check circle2.c -o circle
vampir: circle.c
	#$(VAMPCC) circle2.c -o circle
	vtcc -vt:cc mpicc -std=c99 -ggdb -Wall -Wextra circle2.c -o circle
	mpiexec -np $(NUMPROC) ./circle $(NUMPROC)
	vampir circle.otf
clean:
	rm -f circle
