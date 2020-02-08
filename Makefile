
All: comparaison_multip Statistiques

Statistiques.o:	Statistiques.c
	gcc -Wall -c Statistiques.c -I/usr/local/include/flint -L/usr/local/lib -lgmp -lflint -lm -fopenmp

comparaison_multip.o: comparaison_multip.c
	gcc -Wall -c comparaison_multip.c -I/usr/local/include/flint -L/usr/local/lib -lgmp -lflint -lm -fopenmp

comparaison_multip: comparaison_multip.o
	gcc -Wall -o comparaison_multip comparaison_multip.o -lgmp -lflint -lm -fopenmp

Statistiques: Statistiques.o
	gcc -Wall -o Statistiques Statistiques.o -lgmp -lflint -lm -fopenmp

clear:
	rm *.o
