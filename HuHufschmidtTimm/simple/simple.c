/*
** simple err-demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>

int *
mistake1 ()
{
  int *buf = malloc (sizeof (int) * 6);
  buf[0] = 1;
  buf[1] = 1;
  buf[2] = 2;
  buf[3] = 3;
  buf[4] = 4;
  buf[5] = 5;
  return buf;
}

int *
mistake2 ()
{
  int *buf;
  buf = malloc (sizeof (int) * 4);
  buf[1] = 2;
  return buf;
}

int *
mistake3 ()
{
  // int mistake2_ = 0;  //  unused variable
  int *buf = (int *) mistake2();
  buf[0] = 3;
  return buf;
}

int *
mistake4 ()
{
  int *buf;
  buf = malloc (sizeof (int) * 5);
  buf[0] = 4;
  // free (buf);
  return buf;
}

int
main (void)
{
  /* Modifizieren Sie die folgenden Zeilen nicht */
  int *p[4] = { &mistake1 ()[1], &mistake2 ()[1], mistake3 (), mistake4 () };

  printf ("1 %d\n", *p[0]);
  printf ("2 %d\n", *p[1]);
  printf ("3 %d\n", *p[2]);
  printf ("4 %d\n", *p[3]);

  /* mhh muss hier noch etwas gefreed werden? */
  /* Fügen sie hier die korrekten aufrufe von free() ein */
      /* welcher Pointer war das doch gleich?, TODO: Fixme... :-) */
  free ((void*) p[0] - sizeof(int));  // für mistake1
  free ((void*) p[1] - sizeof(int));  // für mistake2
  free (p[2]);                        // für mistake3
  free (p[3]);                        // für mistake4
  return 0;
}
