#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_vec.h"
#include <math.h>
#include <string.h> 


/*
   construire un polynome en deux variables à partir d'un fichier 
 */
int  read_coeff(char *nom_fichier,fmpz_poly_t * tab){
	
	int i,j;
	int deg1,deg2;
	char buffer[10];
	FILE *F;
	fmpz_t coeff;
	fmpz_init(coeff);
	F=fopen(nom_fichier,"r");
	if(F==NULL){
		printf("erreur d'ouverture");
		return 0;
	}
	fgets( buffer,10,F);	
	sscanf(buffer,"%d \t %d",&deg1,&deg2);	
	for(i=0;i<=deg2;i++){
			fmpz_poly_init2(tab[i],(deg1+1));
			for(j=0;j<=deg1;j++){	
								
				fmpz_fread (F,coeff);
				fmpz_poly_set_coeff_fmpz(tab[i],j,coeff);
				
			
			}
	}
	fmpz_clear(coeff);	
	fclose(F);
	return 0;
}

/*
   génère un polynome modulo un nombre(p) en utilisant la fonction fmpz_mod
 */
void generer_poly_mod(fmpz_poly_t f,fmpz_t p,fmpz_poly_t g){
	fmpz_t coefC;
	fmpz_init(coefC);
	fmpz_t coeff;
	fmpz_init(coeff);
	slong l=fmpz_poly_length (f);
	for(int i=0;i<l;i++){
		fmpz_poly_get_coeff_fmpz(coeff,f,i);
		fmpz_mod(coefC,coeff,p);
		fmpz_poly_set_coeff_fmpz(g,i,coefC);
	}	
	fmpz_clear(coefC) ;
	fmpz_clear (coeff);
	
}
/*
 * renvoie le tableau des coefficients  mod p
 */
 
fmpz_poly_t * tab_poly_mod(fmpz_poly_t * tab, fmpz_t p,int taille){
	fmpz_poly_t  *tab_mod=malloc(sizeof(fmpz_poly_t)*taille+1);
	for(int i=0;i<taille;i++){
		slong l=fmpz_poly_length (tab[i]);
		fmpz_poly_init2(tab_mod[i],(l+1));
		generer_poly_mod(tab[i],p,tab_mod[i]);
	
	}
	if(fmpz_poly_is_zero(tab[taille-1])){
	    printf("le monome du plus haut degre est nul");
	  }
	return tab_mod;
}

	
/*
  *( evaluer un polynome en deux variables (2eme variable) en a )modulo p
  */
void evaluer_1(fmpz_poly_t * f,int taille,fmpz_t premier,fmpz_t a, fmpz_poly_t f_barre){
	fmpz_t res;
	fmpz_init(res);    //on stocke l'évaluation de f[i] en a dans res
	fmpz_t res_mod;
	fmpz_init(res_mod);
	for(int i=0;i<taille;i++){
		fmpz_poly_evaluate_horner_fmpz(res,f[i],a);
		fmpz_mod(res_mod,res,premier);
		fmpz_poly_set_coeff_fmpz(f_barre,i,res_mod);
	}	
}

/*
   evaluer  un polynome en deux variable  en a
 */
void evaluer_2(fmpz_poly_t * f,int taille,fmpz_t a, fmpz_poly_t f_barre){
	fmpz_t res;
	fmpz_init(res);    //on stocke l'évaluation de f[i] en a dans res
	
	for(int i=0;i<taille;i++){
		fmpz_poly_evaluate_horner_fmpz(res,f[i],a);
	
		fmpz_poly_set_coeff_fmpz(f_barre,i,res);
	}	
}
/*  calcul de la borne de hadamard   en utilisant la norme euclidienne
 * 
 */
void borne_H(fmpz_poly_t * f, fmpz_poly_t * g,int taille_f,int taille_g, fmpz_t born){
	slong deg_f;
	slong deg_g;
	fmpz_poly_t norme_f;   //  stocke la norme de chaque f[i]
	fmpz_poly_t norme_g;   //  stocke la norme de chaque g[i]
	fmpz_poly_init(norme_f);
	fmpz_poly_init(norme_g);
	
	fmpz_t NE_f;  // stocke la norme  euclidinne de f
	fmpz_t NE_g;  //stocke la norme euclidienne  de g
	fmpz_init(NE_f);
	fmpz_init(NE_g);
	fmpz_t coeff;
	fmpz_init(coeff);
	/*
	 * calcul de la norme euclidienne de chaque coeff du polynome à deux variables
	 */
	for(int i=0;i<taille_f;i++){	
		fmpz_poly_2norm (coeff,f[i]);
		fmpz_poly_set_coeff_fmpz(norme_f,i,coeff);
	}
	
	for(int i=0;i<taille_g;i++){	
		fmpz_poly_2norm (coeff,g[i]);
		fmpz_poly_set_coeff_fmpz(norme_g,i,coeff);
	}
	
	fmpz_poly_2norm (NE_f,norme_f);
	
	fmpz_poly_2norm (NE_g,norme_g);
	
	deg_f=fmpz_poly_degree(norme_f);
	deg_g=fmpz_poly_degree(norme_g);
	fmpz_pow_ui(NE_f,NE_f,deg_g);
	fmpz_pow_ui(NE_g,NE_g,deg_f);
	fmpz_mul(born,NE_f,NE_g);
	
	
	
	
	// Liberer la mémoire
	
	fmpz_clear(coeff);
	fmpz_clear(NE_f);
	fmpz_clear(NE_g);
	fmpz_clear(coeff);
	fmpz_poly_clear(norme_f);
	fmpz_poly_clear(norme_g);
}

void interpolate(fmpz_poly_t * f, fmpz_poly_t * g, int taille_f,int taille_g, fmpz_t prime){

 
        int borne_res=((taille_f-1)*(taille_g-1))+1;  // borne sur le degré du résultant
	/*
	  Variables temporaires utiles pour le calcul du polynome interpolateur
	 */
	fmpz_poly_t ResPoly;                     //un polynome qui stockera les differentes valeurs du resultant à fin de les caster en un vecteur 
	fmpz_poly_init2(ResPoly,borne_res);
	fmpz_poly_t X_Res;                      //un polynome qui stockera les nombres tirés aléatoirement  afin de les caster en un vecteur 
	fmpz_poly_init2(X_Res,borne_res);

	fmpz_poly_t interpol;                   //polynome interpolateur 
	fmpz_poly_init2(interpol,borne_res+1);
	flint_rand_t state;
	flint_randinit(state);
	fmpz * X= _fmpz_vec_init(borne_res);
	fmpz * Y=_fmpz_vec_init(borne_res);
	fmpz_t * Resultant;                     //tableau qui stockera les differentes valeurs du résultant
	 
	/* on stockera f et g après évaluation en p_evaluate dans des structures nommées f_evaluate et g_evaluate
	 */
	fmpz_poly_t f_evaluate;
	fmpz_poly_t g_evaluate;
	fmpz_poly_init2(f_evaluate,taille_f);
	fmpz_poly_init2(g_evaluate,taille_g);
	
	fmpz_t p_evaluate;  // sert à stocker le nombre aléatoir, qui sera utilisé pour évaluer f et g
	fmpz_init(p_evaluate);	
	
	Resultant=malloc(sizeof(fmpz_t)*(borne_res));
       	fmpz_t * tab=malloc(sizeof(fmpz_t)*(borne_res));  // stockera les a_i tirés aléatoirement (des nombres  tirés aléatoirement pour évaluer f et g)
	
	fmpz_randm(p_evaluate,state,prime);	  
	
	fmpz_init(tab[0]);
	fmpz_set(tab[0],p_evaluate);
	fmpz_poly_set_coeff_fmpz(X_Res,0,p_evaluate);
	_fmpz_vec_set(X,fmpz_poly_get_coeff_ptr(X_Res,0),0);
	//	_fmpz_vec_set(X,fmpz_poly_lead(X_Res),0);
	printf("\n coeff\t ");
	_fmpz_vec_print(fmpz_poly_get_coeff_ptr(X_Res,0),1);

	
	
	evaluer_1(f,taille_f,prime,p_evaluate,f_evaluate);	  //évaluer les deux polynomes en p_evaluate	
	evaluer_1(g,taille_g,prime,p_evaluate,g_evaluate);
	
	fmpz_poly_resultant_euclidean(Resultant[0],f_evaluate,g_evaluate);  //  Calcul du Resultant
	fmpz_mod(Resultant[0],Resultant[0],prime);
	fmpz_poly_set_coeff_fmpz(ResPoly,0,Resultant[0]);
	printf("\n p_evalute: \t");
       	fmpz_print(p_evaluate);
	printf("\t  Resultant: \t");
	fmpz_print(Resultant[0]);
	printf("\n");
	int egal;
	
	for(int i=1;i<borne_res;i++){ 
        
		/*
		 *  verifier que les nombres tirés sont deux à deux distincts
		 */
		
		 do{
			egal=0;
			fmpz_randm(p_evaluate,state,prime);
			for(int j=0;j<i;j++){
				if(fmpz_cmp(p_evaluate,tab[j])==0){
					egal=1;
				}
			}
			
			
		}while(egal==1);		
	       	fmpz_set(tab[i],p_evaluate);
		fmpz_poly_set_coeff_fmpz(X_Res,i,p_evaluate);
		//_fmpz_vec_set(X,fmpz_poly_lead(X_Res),borne_res);
		_fmpz_vec_set(X,fmpz_poly_get_coeff_ptr(X_Res,i),i);
		printf("\ncoeff \t");
		_fmpz_vec_print(fmpz_poly_get_coeff_ptr(X_Res,i),1);
	
		/*
		 * Evaluer les deux polynomes f et g en p_evaluate
		 */
		evaluer_1(f,taille_f,prime,p_evaluate,f_evaluate);		
		evaluer_1(g,taille_g,prime,p_evaluate,g_evaluate);
		/*
		 * je calcule le resultant
		 */
		fmpz_init(Resultant[i]);
	        fmpz_poly_resultant_euclidean(Resultant[i],f_evaluate,g_evaluate);
	       	fmpz_mod(Resultant[i],Resultant[i],prime);
		fmpz_poly_set_coeff_fmpz(ResPoly,i,Resultant[i]);
	        printf("\n p_evalute: \t");
		fmpz_print(p_evaluate);
		printf("\t  Resultant: \t");
		fmpz_print(Resultant[i]);
		
	}
	printf("\n %d\n",borne_res);
	fmpz_print(prime);
	printf("\n vexteurX: ");
       	fmpz_poly_print(X_Res);
	//_fmpz_vec_set(X,fmpz_poly_lead(X_Res),borne_res);
	//	X=fmpz_poly_lead(X_Res);
	//	Y=fmpz_poly_lead(ResPoly);
       	printf("\n vecteur \t");
       	_fmpz_vec_print(X,borne_res);
	//fmpz_poly_interpolate_fmpz_vec(interpol,X,Y,borne_res);
	//	fmpz_poly_print(interpol);


	/*
	 * liberation de  la mémoire 
	 */
	for(int k=0;k<borne_res;k++){
	  fmpz_clear(Resultant[k]);
	  fmpz_clear(tab[k]);
	}
	fmpz_clear(p_evaluate);
	flint_randclear(state);
	free(tab);
	free(Resultant);
	fmpz_poly_clear(f_evaluate);
	fmpz_poly_clear(g_evaluate);
	fmpz_poly_clear(interpol);
	fmpz_poly_clear(ResPoly);
	fmpz_poly_clear(X_Res);
	_fmpz_vec_clear(X,borne_res);
	_fmpz_vec_clear(Y ,borne_res);
	
} 
fmpz_t *  generer_nb_premier(fmpz_t Born,int bits_premier,int * nb){
  flint_rand_t state;
  flint_randinit(state);
  slong r;
  fmpz_t mi;
  fmpz_init(mi);
  fmpz_t M;
  fmpz_one(M); //stockera le produit des nombres premiers

   /*
   * calcul du produit des nombres premiers 
   */
  fmpz_mul_ui(Born,Born,2);
  fmpz_add_ui (Born,Born,1);
  r=fmpz_clog_ui(Born, 2);
  *nb=r;
  fmpz_t * tab;   //stockera les nombres premiers 
  int egal;
  do{	
    
    fmpz_t tmp;
    fmpz_init(tmp);
    tab=malloc(sizeof(fmpz_t)*r);
    for(int i=0;i<r;i++){
      do{
	egal=0;
	fmpz_randbits(mi,state,bits_premier);
	if(i>0){
	  for(int j=0;j<i;j++){
	    if(fmpz_cmp(mi,tab[j])==0){
	      egal=1;
	    }
	  }
	}
      }while(fmpz_is_probabprime_lucas(mi)==0 || egal==1);
      fmpz_mul(M,M,mi);
      fmpz_init_set(tab[i],mi);		
      
    }
    
 
  }while((fmpz_cmp(M,Born)<0));//|| (fmpz_cmp(M,dM)>0

   /*
    libérer la mémoire
  */
  flint_randclear(state);
  fmpz_clear(M);
  fmpz_clear(mi);
  return tab; 
}

void restes_chinois(fmpz_poly_t * f,fmpz_poly_t * g,int taille_f,int taille_g,int bits_premier,int *nb){
  
  /*
    Recuperer le nombre exacte des nombres premiers
   */
  fmpz_t Born;  //borne de hadamard
  int r;// le nombre exacte des nombres premiers
  fmpz_t * tab;
  borne_H(f,g,taille_f,taille_g,Born);
  tab=generer_nb_premier( Born, bits_premier,&r);

  /*
    des variables temporaires utiles pour le théorème des restes chinois
   */
  fmpz_t r1;  
  fmpz_t r2; 
  fmpz_t m1;
  fmpz_t m2;
  fmpz_t x; //pour stocker la solution
  fmpz_init(m1);
  fmpz_init(m2);

  
  flint_rand_t state;
  flint_randinit(state);
 
  fmpz_t resultant;
  fmpz_init(resultant);
  /*
   * évaluer f et g en un nombre tiré aléatoirement et les stocker respectivement  dans des polynomes à une seule variable  f_evaluate et g_evaluate
   
   */
  fmpz_poly_t f_evaluate;  
  fmpz_poly_t g_evaluate;
  fmpz_poly_init2(f_evaluate,taille_f);
  fmpz_poly_init2(g_evaluate,taille_g);

  
  /*
   *  f mod p_i et g mod p_i  ( p_i c'est les nombres premiers calculés précédemment)
    vérifier que le monome du plus haut degé n'est pas nul 
   */
  for(int k=0;k<r;k++){
    fmpz_poly_t * f_mod=tab_poly_mod(f,tab[k],taille_f);
    fmpz_poly_t * g_mod=tab_poly_mod(g,tab[k],taille_g);
  }
  
  /*
   *  appliquer le théorème des restes chiniois 
   */
	
  /*
   * générer un nombre aléatoire pour evaluer f et g
   */
 
  fmpz_t eval;
  fmpz_randbits(eval,state,bits_premier);
  
  evaluer_2(f,taille_f,eval,f_evaluate);
  evaluer_2(g,taille_g,eval,g_evaluate);
  /*
    calcul du resultatnt modulo m0 et appliquer le théorème des restes chinois
   */
  fmpz_t mod_res; // pour stocker resultant mod M(produit des nombres premiers )
  fmpz_poly_resultant_euclidean(resultant,f_evaluate,g_evaluate); // calcul du resultant des deux polynomes
  printf("resultant\t");
  fmpz_print(resultant);
  printf("\n");
 
  fmpz_mod(r1,resultant,tab[0]);
  
  fmpz_init_set(m1,tab[0]);
  for(int l=1;l<r;l++){
     
      fmpz_mod(r2,resultant,tab[l]);    
      fmpz_init_set(m2,tab[l]);
      fmpz_CRT(x,r1,m1,r2,m2,0);
      fmpz_mul(m1,m1,m2);     
      fmpz_init_set(r1,x);
      fmpz_mod(r1,r1,m1);
      
    
  }
  fmpz_print(r1);
  printf("\n");
  fmpz_mod(mod_res,resultant,m1);
  fmpz_print(mod_res);
 
	
  /*
   * Liberation de la memoire
   */
  flint_randclear(state);
  fmpz_clear(m1);
  fmpz_clear(m2);
  fmpz_clear(r1);
  fmpz_clear(r2);
  fmpz_clear(x);
  fmpz_clear(mod_res);
  free(tab);
  fmpz_poly_clear(f_evaluate);
  fmpz_poly_clear(g_evaluate);
	
	
}
 

 
int main(int argc,char *argv[]){
	
	/*
	 * Récupérer la borne de hadamard   et la stocker dans born afin de calculer le cardinal des nombres premiers à tirer aléatoirement 
	 */
	 fmpz_t born;	 
	 fmpz_init(born);
	 int nombre_p; 
	
	 /*
	   récupérer le degré des polynomes à partir d'un fichier 
	  */
	char buff[10];
	int deg1_f,deg2_f;
	int deg1_g,deg2_g;
	int bits_premier=atoi(argv[1]); //je récupère la taille du nombre premier
	flint_rand_t state;
	flint_randinit(state);
	fmpz_t p;  // nombre premier 
	/*
	 * générer le nombre premier p
	 */
	fmpz_init(p);
	do{
		fmpz_randbits(p,state,bits_premier);
	}while( fmpz_is_probabprime_lucas(p)==0);
	
	
	/*
	 * Lire les coeff des deuc polynomes à deux variables f et g à partir des fichiers F et G
	 */
	FILE * F;
	FILE * G;
	G=fopen("fichier_2.txt","r");
	F=fopen("new_fichier.txt","r");
	if(G==NULL){
		printf("  erreur d'ouverture de G");
		return 0;
	}
	if(F==NULL){
		printf(" erreur d'ouverture de F");
		return 0;
	}
	/*
	 * générer le  polynome à deux variables f
	 */
	fgets(buff,10,F);	
	sscanf(buff,"%d \t %d",&deg1_f,&deg2_f);
	fmpz_poly_t * f;
	f=malloc(sizeof(fmpz_poly_t)*(deg2_f+1));      //tableau qui stocke le polynome à deux variables 
	read_coeff("new_fichier.txt",f);			//je construit le polynome à deux variable à partir d'un fichier 
	/*
	 * générer le  polynome à deux variables g
	 */
	fgets(buff,10,G);
	fmpz_poly_t * g;
	sscanf(buff,"%d \t %d",&deg1_g,&deg2_g);
	g=malloc(sizeof(fmpz_poly_t)*(deg2_g+1));     //tableau qui stocke le polynome à deux variables 
	read_coeff("fichier_2.txt",g);			   //je construit le polynome à deux variable à partir d'un fichier 
	
	/*
	 *  Récuperer la borne d'Hadamard
	 */
	 borne_H(f,g,deg2_f+1,deg2_g+1,born);
	 printf(" \n borne de hadamard est \n");
	 fmpz_print(born);
	 printf("\n");
	/*
	 * Methode du prof pour calculer le produit des nombres premiers 
	 *  fmpz_t nb_premier;
	 fmpz_fdiv_q_si(nb_premier,born,bits_premier);
	 fmpz_add_ui(nb_premier,nb_premier,1);
	 slong nb=fmpz_get_si (nb_premier);
	 flint_printf("\n nb_premier: %wd \t  \n",nb);
	 
	 
	 /*
	  *  probleme d'interpolation 
	  */
      	 interpolate(f,g,deg2_f+1,deg2_g+1,p);
	
	/*
	 * Theoreme des restes chinois
	 */
	// restes_chinois(f,g,deg2_f+1,deg2_g,bits_premier,&nombre_p);
		// generer_nb_premier(born,bits_premier);
	

	 // Libérer la mémoire dynamique

	for(int i=0;i<=deg2_f;i++){
	  fmpz_poly_clear(f[i]);	
	  fmpz_poly_clear(g[i]);	
	}
	fmpz_clear(born);
	free(f);
	free(g);
	fmpz_clear(p);
	flint_randclear(state);
	fclose(F);
	fclose(G);
	return 0;
}
					
