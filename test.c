#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include <ulong_extras.h>
#include <math.h>
#include <string.h>


int creer_fichier(char *nom_fichier,int n1,int n2){
	FILE *f;
	fmpz_t x;
	fmpz_init(x);
	flint_rand_t state;
	flint_randinit(state);
	f=fopen(nom_fichier,"w");
	if(f==NULL){
		printf("erreur d'ouverture");
		return 0;
	}	
	for(int i=0;i<=n2;i++){
		for(int j=0;j<=n1;j++){
			fmpz_randtest(x,state,128);
			fmpz_fprint(f,x);
			fprintf(f,"\t");
		}
		fprintf(f,"\n");
			
	}
	printf("\n");	
	fmpz_clear(x);
	flint_randclear(state);
	return (n1+1)*(n2+1);
}
int main(){
	creer_fichier("new_fichier.txt",10,12);
	return 0;
}
