/**
 * @file HyperTan_custom.cpp
 * HYPERTAN activation function
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>, Eleonora.
 * @version 0.1
 */

#include "HyperTan_custom.h"

 HyperTan_custom::HyperTan_custom(Bundle *B, int dim_sys, double* coeff, int deg){


   // Initialize model
   strcpy(this->name,"HyperTan3D");

   // List of state variables
   symbol x1("x1"), x2("x2"), x3("x3");
   lst vars;
   lst dyns;

   // System's dynamics
   ex dx1 = 0;
   ex dx2 = 0;
   ex dx3 = 0;
   
   for(int i=0; i<=deg; i++){
   	dx1 = dx1 + coeff[i]*pow(x1,i);
	dx2 = dx2 + coeff[i]*pow(x2,i);
	dx3 = dx3 + coeff[i]*pow(x3,i);
   }

   if(dim_sys==1){
   	vars = x1;
	dyns = dx1;
   } else if(dim_sys==2){
	vars = {x1, x2};
	dyns = {dx1,dx2};
   } else if(dim_sys==3){
	vars = {x1, x2, x3};
	dyns = {dx1,dx2,dx3};
   }

   this->vars = vars;
   this->dyns = dyns;


   // Init reach set
   this->reachSet = B;

}
