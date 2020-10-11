/**
 * @file HyperTan3D.cpp
 * HYPERTAN activation function
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>, Eleonora.
 * @version 0.1
 */

#include "HyperTan3D.h"

 HyperTan3D::HyperTan3D(Bundle *B, int dim_sys){


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
   
   //minmax approx with deg=3
   vector< double > coeff {0.9266, -0.1854, 0.0239, -0.001494, 0.00003511};
   dx1 = dx1 + coeff[0]*pow(x1,1) + coeff[1]*pow(x1,3) + coeff[2]*pow(x1,5) + coeff[3]*pow(x1,7) + coeff[4]*pow(x1,9);
   dx2 = dx2 + coeff[0]*pow(x2,1) + coeff[1]*pow(x2,3) + coeff[2]*pow(x2,5) + coeff[3]*pow(x2,7) + coeff[4]*pow(x2,9);
   dx3 = dx3 + coeff[0]*pow(x3,1) + coeff[1]*pow(x3,3) + coeff[2]*pow(x3,5) + coeff[3]*pow(x3,7) + coeff[4]*pow(x3,9);


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
