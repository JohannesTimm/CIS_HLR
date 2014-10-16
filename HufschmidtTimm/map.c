#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Definieren Sie ein enum cardd
typedef enum { X = 0, N = 1, S = 2, E = 4, W = 8} cardd;

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
cardd* map[3][3] = {{X, X, X}, {X, X, X}, {X, X, X}};

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir)
{
  int d = (int) dir;
  if ((x < 0) | (x > 2)) {
    printf ("Ungültiger Wert %d für x\n", x);
    // exit(EXIT_FAILURE);
    return;
  }
  if ((y < 0) | (y > 2)) {
    printf ("Ungültiger Wert %d für y\n", y);
    // exit(EXIT_FAILURE);
    return;
  }
  if ((d < 0) | (d > 15) | (d == 3) | (d == 12)) {
    printf ("Ungültiger Wert %d für d\n", d);
    // exit(EXIT_FAILURE);
    return;
  }
  map[x][y] = dir; // Kompilierungs-Fehler
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
  int row, col;
  cardd dir;
  char* singleLetters = "NSEW";
  for (row = 0; row < 2; row++) {
    for (col = 0; col < 2; col++) {
      dir = map[row][col]; // Kompilierungs-Fehler
      switch (dir) {
        case N:
          //
        break;
        case S:
        //
        break;
        case E:
        //
        break;
        case W:
        //
        break;
        default:
        //
        break;
      }
    }
  }
}

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
