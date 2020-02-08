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
#include <math.h>

int main(int argc, char ** argv){
    if(argc!=3){
       fprintf(stderr, "Bad usage !!\n Use: %s <deg> <nb_bits_coeff_max>\n", argv[0]);
       return 1;
    }
    unsigned long int deg = atoi(argv[1]);
    printf ("le degres du polynome est: %ld\n", deg);

    unsigned long int nb_bits_coeff_max= atoi(argv[2]);

    flint_rand_t state;
    flint_randinit(state);
	// déclaration
	fmpz_poly_t x, y, z;
	// initialisation ( allocation mémoire )
	fmpz_poly_init2(x,deg);
	fmpz_poly_init2(y,deg);
	fmpz_poly_init2(z,deg*2+1);



	//génération aléatoire
	fmpz_poly_randtest(x,state, deg,nb_bits_coeff_max);
	fmpz_poly_randtest(y, state, deg, nb_bits_coeff_max);
	
	
	//Tests de comparaison
	

	double total1=0;
	double total2=0;
      	double total3=0;
        double total4=0;
       	double total6=0;
       	double total7=0;
	double total8=0;
       	double total9=0;
	int j;
	double e,e1;
	
	for(j=0; j<15000; j++){

	  e=omp_get_wtime();
	  fmpz_poly_mul(z, x, y);
	  e1=omp_get_wtime()-e;
	  printf( "e= %0.12lf \n", e);
	  fprintf(stdout, "time mul:%0.12lf secondes \n", e1);
	  total1+=e1;
	}
	fprintf(stdout, "time mul total:%12lf secondes\n", total1);
	
	//Résultats des tests:
	/*

	  Degré 500, nb_bits_max= 50
	  time mul total:    3.242704 secondes

	   Degré 500, nb_bits_max= 56
	  Test 4: time mul total:    3.497484 secondestime 

	  Degré 500, nb_bits_max: 120
	  time mul total:    2.725789 secondes

	  Degré: 500, nb_bits_max: 200
	  time mul total:   14.871770 secondes

	  Degré 500 nb_bits_max: 300
	  time mul total:   20.305168 secondes
	  
	  Degré 1000, nb_bits_max: 56
	  Test 3: time mul total:    8.735670 secondes

	   Degré 1000 nb_bits_max=70
	  time mul total:   10.946877 secondes

	  Degré 1001 nb_bits_max 80
	  time mul total:   12.580017 secondes

	  Degré 2000 nb_bits_max: 90
	  time mul total:   37.919530 secondes*/


	 for(j=0; j<15000; j++){

	    e=omp_get_wtime();
	    fmpz_poly_mul_classical(z, x, y);
	    e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mul_classical:%lf secondes \n", e1);
	    total2+=e1;
	  }
	  fprintf(stdout, "mul_classical total:%lf secondes", total2);
	
	  /*
	    degré 500 nb_bits_max: 50
	    mul_classical total:6.089227 secondes

	    degré 500 nb_bits_max: 200
	    mul_classical total:12.263215 secondes
	    
	    degré: 500 nb_bits_max=300
	    mul_classical total:14.710557
	    mul_classical total:14.756642 secondes

	    degré 500 nb_bits 56
	    mul_classical total:6.349693 secondes
	    
	    degré 1000 nb_bits_max: 70
	    mul_classical total:28.917602 secondes

	    degré 1000 nb_bits_max=56
	    mul_classical total:25.692584 secondes
	    
	    Degré 1001 nb_bits_max 80
	    mul_classical total:31.206057 secondes

	     Degré 2000 nb_bits_max: 90
	     mul_classical total:139.373415 secondes

ygbyg

	  */
	  
	  for(j=0; j<15000; j++){

	  e=omp_get_wtime();
	  fmpz_poly_mullow_classical(z, x, y,fmpz_poly_length(x)+fmpz_poly_length(y)-1 );
	  e1=omp_get_wtime()-e;
	  fprintf(stdout, "time mullow_classical:%lf secondes \n", e1);
	  total3+=e1;
	  }
	  fprintf(stdout, "mullow_classical total:%lf secondes", total3);
	  
	  
	  for(j=0; j<15000; j++){

	    e=omp_get_wtime();
	    fmpz_poly_mulhigh_classical(z, x, y, 0);
	    e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mulhigh_classical:%lf secondes \n", e1);
	    total4+=e1;
	  }
	  fprintf(stdout, "mulhigh_classical total:%lf secondes", total4);
	  
	  for(j=0; j<15000; j++){

	    e=omp_get_wtime();
	    fmpz_poly_mul_karatsuba(z, x, y);
	    e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mul_karatsuba:%lf secondes \n", e1);
	    total6+=e1;
	  }
	  fprintf(stdout, "time mul_karatsuba total:%lf secondes", total6);
	  
	  for(j=0; j<15000; j++){

	    e=omp_get_wtime();
	    fmpz_poly_mullow_karatsuba_n(z, x, y, fmpz_poly_length(x)+fmpz_poly_length(y)-1);
	    e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mullow_karatsuba_n:%lf secondes \n", e1);
	    total7+=e1;
	  }
	  fprintf(stdout, "time mullow_karatsuba_n total:%lf secondes", total7);
	  
	  
	   for(j=0; j<15000; j++){

	     e=omp_get_wtime();
	    fmpz_poly_mullow(z,x,y, fmpz_poly_length(x)+fmpz_poly_length(y)-1);
	    e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mullow :%lf secondes \n", e1);
	    total8+=e1;
	   }
	   fprintf(stdout, "time mullow total:%lf secondes", total8);
   

   
	   for(j=0; j<15000; j++){

	     e=omp_get_wtime();
	    fmpz_poly_mulhigh_n(z, x, y, fmpz_poly_length(x)+fmpz_poly_length(y)-1);
	    e1=omp_get_wtime()-e;
	    fprintf(stdout, "time mulhigh:%lf secondes \n", e1);
	    total9+=e1;
	   }
	   fprintf(stdout, "time mulhigh total:%lf secondes\n", total9);
	   
	   fprintf(stdout, "mul: %lf \n", total1);
	   fprintf(stdout, "mul_classical: %lf \n", total2);
	   fprintf(stdout, "mullow_classical: %lf \n", total3);
	   fprintf(stdout, "mulhigh_classical: %lf\n", total4);
	   fprintf(stdout, "mul_karatsuba: %lf\n", total6);
	   fprintf(stdout, "mullow_karatsuba_n: %lf \n", total7);
	   fprintf(stdout, "mullow: %lf\n", total8);
	   fprintf(stdout, "mulhigh: %lf \n", total9);
}
