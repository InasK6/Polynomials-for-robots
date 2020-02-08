// comparaison entre les differentes fonctions de multiplication  de FLINT
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include <math.h>
#include <omp.h>
#include "ulong_extras.h"
#include <string.h>
int main(int argc, char ** argv){
    if(argc!=3){
       fprintf(stderr, "Bad usage !!\n Use: %s <deg> <coeff_max>\n", argv[0]);
       return 1;
    }
    unsigned long int deg = atoi(argv[1]);
    printf ("le degres du polynome est: %ld\n", deg);

    unsigned long int coeff_max= atoi(argv[2]);
    int i=0;
    while(coeff_max>0){
      coeff_max/=2;
      i++;
    }

    flint_rand_t state;
    flint_randinit(state);
	// déclaration
	fmpz_poly_t x, y, z;
	// initialisation ( allocation mémoire )
	fmpz_poly_init2(x,deg);
	fmpz_poly_init2(y,deg);
	fmpz_poly_init2(z,deg*2+1);



	//génération aléatoire
	fmpz_poly_randtest(x,state, deg,i);

	fmpz_poly_randtest(y, state, deg, i);

    /*
    * génère les mêmes polynômes après chaque execution est-ce que c'est normal?
    * ne marche pas à 10000 ( n'alloue pas le deuxième polynôme
    * va jusqu'à 5000
    */  	
	
	//Test 	
	

	double e, e1;
    double total1=0;
    double total2=0;
    double total3=0;
    double total4=0;
    double total5=0;
    double total6=0;
    double total7=0;
    double total8=0;
    double total9=0;
	// multiplication
	int j;
	//1)
   for(j=0; j<1000; j++){

        e=omp_get_wtime();
	    fmpz_poly_mul(z, x, y);
        e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mul:%lf secondes \n", e1);
        total1+=e1;
    }
    fprintf(stdout, "time mul total:%lf secondes", total1);
    /* pour 10000 : 8 secondes
       pour 1000: 0,9 secondes
    */

	// 2)
   /* for(j=0; j<1000; j++){

        e=omp_get_wtime();
	    fmpz_poly_mul_classical(z, x, y);
        e1=omp_get_wtime()-e;
        fprintf(stdout, "time mul_classical:%lf secondes \n", e1);
        total2+=e1;
    }
    fprintf(stdout, "mul_classical total:%lf secondes", total2);*/
    /*
    mul_classical total:240.622548 secondes pour 1000 itération
    */

		
	//3)
	// deg(x*y)=deg(x)+deg(y) et length=deg+1
     /*for(j=0; j<1000; j++){

        e=omp_get_wtime();
	    fmpz_poly_mullow_classical(z, x, y,fmpz_poly_length(x)+fmpz_poly_length(y)-1 );
        e1=omp_get_wtime()-e;
        fprintf(stdout, "time mullow_classical:%lf secondes \n", e1);
        total3+=e1;
    }
    fprintf(stdout, "mullow_classical total:%lf secondes", total3);*/
    /* pour 1000 itérations: mullow_classical total:241.787782 secondes
    */

	//4)
     /*for(j=0; j<1000; j++){

        e=omp_get_wtime();
	    fmpz_poly_mulhigh_classical(z, x, y, 0);
        e1=omp_get_wtime()-e;
        fprintf(stdout, "time mulhigh_classical:%lf secondes \n", e1);
        total4+=e1;
    }
   fprintf(stdout, "mulhigh_classical total:%lf secondes", total4);*/
    /* mulhigh_classical total:238.771351 secondes pour 1000 itérations
    */

	//5) -> -> Problème avec le résultat <- <-
   /* e=omp_get_wtime();
	fmpz_poly_mulmid_classical(z, x, y);*/
	// ATTENTION: assumes len1>=len2
	//fprintf(stdout, "time fmpz_poly_mulmid:%lf secondes \n", omp_get_wtime()-e);
	// Resultat: time fmpz_poly_mulmid:0.000016 secondes 
	/*Sets res to the middle len(poly1)- len(poly2)+ 1 coefficients of poly1 * poly2, i.e.
	the coefficient from degree len2 - 1 to len1 - 1 inclusive. Assumes that len1 >=
	len2.*/
	//Elle renvoie 0, n'effectue pas bien la multiplication!! je n'ai pas compris sa définition
	
	//6) 
    /*for(j=0; j<1000; j++){

        e=omp_get_wtime();
	    fmpz_poly_mul_karatsuba(z, x, y);
        e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mul_karatsuba:%lf secondes \n", e1);
        total6+=e1;
    }
    fprintf(stdout, "time mul_karatsuba total:%lf secondes", total6);*/
    /*
    * time mul_karatsuba total:100.985326 secondes
    */

	//7)
    /*for(j=0; j<1000; j++){

        e=omp_get_wtime();
	    fmpz_poly_mullow_karatsuba_n(z, x, y, fmpz_poly_length(x)+fmpz_poly_length(y)-1);
        e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mullow_karatsuba_n:%lf secondes \n", e1);
        total7+=e1;
    }
    fprintf(stdout, "time mullow_karatsuba_n total:%lf secondes", total7);*/
     /*time mullow_karatsuba_n total: 285.093917 
    */


	//8)
    for(j=0; j<15000; j++){

        e=omp_get_wtime();
	    fmpz_poly_mullow(z,x,y, fmpz_poly_length(x)+fmpz_poly_length(y)-1);
        e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mullow :%lf secondes \n", e1);
        total8+=e1;
    }
    fprintf(stdout, "time mullow total:%lf secondes", total8);
    /*time mullow total:0.894722 secondes
    */

	//9) --> demander si c'est bien ce que j'ai compris pour n <--
    for(j=0; j<15000; j++){

        e=omp_get_wtime();
	    fmpz_poly_mulhigh_n(z, x, y, fmpz_poly_length(x)+fmpz_poly_length(y)-1);
        e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mulhigh:%lf secondes \n", e1);
        total9+=e1;
    }
    fprintf(stdout, "time mulhigh total:%lf secondes", total9);
    /*
    * time mulhigh total:0.921023 secondes
    */

    
	
    
	 

	
	
	//affichage
    printf("polynôme 1: ");
	fmpz_poly_print(x); printf("\n");
     printf("polynôme 2: ");
	fmpz_poly_print(y); printf("\n");
    printf("résultat de multiplication: ");
	fmpz_poly_print(z); printf("\n");

    /*
    Test 2 avec 1000 itérations
    */
    fprintf(stdout, "mul: %lf \n", total1);
     fprintf(stdout, "mul_classical: %lf \n", total2);
     fprintf(stdout, "mullow_classical: %lf \n", total3);
     fprintf(stdout, "mulhigh_classical: %lf\n", total4);
     fprintf(stdout, "mul_karatsuba: %lf\n", total6);
     fprintf(stdout, "mullow_karatsuba_n: %lf \n", total7);
     fprintf(stdout, "mullow: %lf\n", total8);
     fprintf(stdout, "mulhigh: %lf \n", total9);
    /* résultat:
    * mul: 0.911379
    *mul_classical: 241.203430
    *mullow_classical: 245.381320
    * mulhigh_classical: 240.954911
    * mul_karatsuba: 98.477235
    * mullow_karatsuba: 282.954468 
    *mullow: 0.853415
    * mulhigh: 0.872971
    */	
    
    /* On refait les tests sur mul, mullow karatsuba, mullow et mulhigh avec plus d'itérations pour déterminer la plus rapide d'entre elles 
    Résultat: 
    
    test 1: 10000 itérations
    mul: 9.186308 
    mullow: 9.472586
    mulhigh: 9.480566 

    test 2: 10000 itérations
    mul: 10.308515 
    mullow: 9.351091
    mulhigh: 9.590245 
    
    test 3: 15000 itérations
    mul: 13.687809
    mullow: 13.365484
    mulhigh: 13.381107 

    test 4: 15000 itérations
   
    mul: 13.558455 
    mullow: 13.321632
    mulhigh: 13.388296

    test 5=15000 itération
    mul: 13.659625 
    mullow: 13.516147
    mulhigh: 13.619575 






*/

    

	//libération de la mémoire
    
	fmpz_poly_clear(x);
   
	fmpz_poly_clear(y);
   
	fmpz_poly_clear(z);

	
	
	return 0;
}


 
// fmpz_poly_t is nomalised if length=0 or the leading coefficient is non-zero 
// ui : unsigned integer
// si: signed integer 
