#include <stdio.h>
#include <stdint.h>

//générer un polynôme modulo un nombre premier type fmpz_mod_poly
void generer_poly_mod(fmpz_poly_t f,fmpz_t p,fmpz_mod_poly g){
    fmpz_t coeff_simple;
    slong n;
    slong length=fmpz_poly_length ( f );
    for( n=0; n<length ;n++){
        fmpz_poly_get_coeff_fmpz(coeff_simple, f, n);
    }

}

// fonction qui génère un nombre premier



int main(){
//polynômes à deux variables: tableau de polynômes à deux variables (x,y)
  
  int d;
  // variable qui définit le degré du polynôme à une variable
  printf("Entrez le degrès du polynôme par rapport à la première variable y");
  scanf(" %d", d);
  // On génére un tableau de polynôme à une variable x de taille d
  fmpz_poly_t poly1[d];
  fmpz_poly_t poly2[d];
  
  
   int i;
   flint_rand_t state;
   flint_randinit(state)
     int bits; // nombre de bits maximal des coefficients
   int nb_bits=0;
    printf("Veuillez rentrer la valeur maximal possible de vos coefficients %d", bits);
    scanf(" %d", bits);
    while(bits>0){
      bits/=2;
      nb_bits++;
    }

    /*
      La stratégie multi-modulaire pour le calcul du résultant s'applique sur des matrices carrées c'est pour ça que je prends le même degrès pour les polynômes à une variable aussi 
     */
  for (i=0; i<d; i++){
    fmpz_poly_init2(poly1[i]);
    fmpz_poly_init2(poly2[i]);
    fmpz_poly_randtest(poly1[i],state, d, nb_bits);
    fmpz_poly_randtest(poly2[i],state, d, nb_bits); 
    
  }

  // calcul de la valeur maximale des valeurs absolues des coefficients des polynômes
  
  // Je récupère la valeur de chaque coefficient et je trouve le max
  fmpz_t B=0;
  
  for (i=0; i<d; i++){
    int j;
    fmpz_t x;
    for(j=0; j<d; j++){
      fmpz_poly_get_coeff_fmpz(x, poly1[i], d);
      if(x>B){
	B=x;
      }
    }
  }
  // calcul de la borne d'Hadamard
  double H= n**(n/2) * B**n

  /*On cherche à calculer des nombres premiers tels que leur produit soient plus grand que la borne d'Hadamard.
   condition: 2⁵⁰< nb_premier<2⁶²
   Je génère d'abord des un fmpz aléatoirement entre les valeurs voulues , je teste s'il est premier sinon j'en génère un autre 
    int fmpz_is_prime ( const fmpz_t n )
    /*
    * The function assumes that n is likely prime, i.e. it is not very efficient if n is composite.
    * A strong probable prime test should be run first to ensure that n is probably pri
    */
   */


    
    // évaluation des polynômes sur les nombres premiers pour la première variable

    // calcul du résultant 
     
  return 0;
}
