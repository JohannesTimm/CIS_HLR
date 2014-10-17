#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

// Definieren Sie ein enum cardd
typedef enum { X = 0, N = 1, S = 2, E = 4, W = 8} cardd;

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
cardd map[3][3] = {{X, X, X}, {X, X, X}, {X, X, X}};

bool valid_dir(cardd dir) {
  bool ok;
  ok = (dir == N) || (dir == S) || (dir == E) || (dir == W) ||
       (dir == (N | W)) || (dir == (S | W)) ||
       (dir == (N | E)) || (dir == (S | E));
  return ok;
}

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir)
{
  int d = (int) dir;
  if ((x < 0) | (x > 2)) {
    printf ("Ungültiger Wert %d für x\n", x);
    return;
  }
  if ((y < 0) | (y > 2)) {
    printf ("Ungültiger Wert %d für y\n", y);
    return;
  }
  if (!valid_dir(dir)) {
    printf ("Ungültiger Wert %d für d\n", d);
    return;
  }
  // printf("x = %d, y = %d, (int)dir = %d\n", x, y, d); // only for test
  map[x][y] = dir;
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void) {
  int row, col, pos, num;
  char const_zeile[] = "         ", zeile[10];
  for (row = 0; row < 3; row++) {
    strcpy(zeile, const_zeile);  // Leere Zeile
    pos = 0;
    for (col = 0; col < 3; col++) {
      num = (int) map[row][col];
      switch (num) {
        case 0:
          zeile[pos] = 'O';
          break;
        case 1:
          zeile[pos] = 'N';
          break;
        case 2:
          zeile[pos] = 'S';
          break;
        case 4:
          zeile[pos] = 'E';
          break;
        case 8:
          zeile[pos] = 'W';
          break;
        case 5:
          if (pos == 0) {
            zeile[pos] = 'N';
            zeile[pos + 1] = 'E';
          } else {
            zeile[pos-1] = 'N';
            zeile[pos] = 'E';
          }
          break;
        case 9:
          if (pos == 0) {
            zeile[pos] = 'N';
            zeile[pos + 1] = 'W';
          } else {
            zeile[pos-1] = 'N';
            zeile[pos] = 'W';
          }
          break;
        case 6:
          if (pos == 0) {
            zeile[pos] = 'S';
            zeile[pos + 1] = 'E';
          } else {
            zeile[pos-1] = 'S';
            zeile[pos] = 'E';
          }
          break;
        case 10:
          if (pos == 0) {
            zeile[pos] = 'S';
            zeile[pos + 1] = 'W';
          } else {
            zeile[pos-1] = 'S';
            zeile[pos] = 'W';
          }
          break;
      } //  end switch
      /* Just for testing
      printf ("row = %d, col = %d, d = %d, dir = ", row, col, num);
      if (pos == 0) {
        printf("%c%c\n", zeile[pos], zeile[pos+1]);
      } else {
        printf("%c%c\n", zeile[pos-1], zeile[pos]);
      }
      */
      pos = pos + 4;
    }  // end for col
    printf("%s\n", zeile);
  }  // end for row
}  // end show_map

int main (void)
{
	// In dieser Funktion darf nichts verändert werden!
	set_dir(0, 1, N);
	set_dir(1, 0, W);
	set_dir(1, 4, W);
	set_dir(1, 2, E);
	set_dir(2, 1, S);

	set_dir(0, 0, N|W);
	set_dir(0, 2, N|E);
	set_dir(0, 2, N|S);
	set_dir(2, 0, S|W);
	set_dir(2, 2, S|E);
	set_dir(2, 2, E|W);

	show_map();

	return 0;
}
