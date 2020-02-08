#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_poly_mat.h"
#include <math.h>
#include <string.h> 

int creer_fichier(char *nom_fichier,int n1,int n2,int bits){
	FILE *f;
	fmpz_t x;
	fmpz_init(x);
	flint_rand_t state;
	flint_randinit(state);
	f=fopen(nom_fichier,"w");
	if(f==NULL){
		printf("erreur d'ouverture");
		return EXIT_FAILURE ;
	}
	fprintf(f,"%d\t %d\n",n1,n2); 
	for(int i=0;i<=n2;i++){
		for(int j=0;j<=n1;j++){
			
			fmpz_randtest_not_zero(x,state,bits);
			fmpz_fprint(f,x);
			fprintf(f,"\t");
		}
		fprintf(f,"\n");
			
	}
	printf("\n");	
	fmpz_clear(x);
	flint_randclear(state);
	fclose(f);
	return 0; 
}
int main(int argc, char *argv[]){
	int deg1=atoi(argv[1]);// degre du polynome à une variable
	int deg2=atoi(argv[2]); // degre du polynome à deux variables
	int bits=atoi(argv[3]);  // la taille binaire des coeff
	creer_fichier("fichier_2.txt",deg1,deg2,bits);
	creer_fichier("new_fichier.txt",deg1,deg2,bits+1);
	return 0;
}
