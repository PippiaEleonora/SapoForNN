/**
 * @file HyperTan.cpp
 * HYPERTAN activation function
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>, Eleonora.
 * @version 0.1
 */

#include "HyperTan.h"

 HyperTan::HyperTan(bool cheb){


   // Initialize model
   strcpy(this->name,"HyperTan");
   int dim_sys = 2;//3;
   // List of state variables
   symbol x1("x1"), x2("x2");//, x3("x3");
   lst vars;
   vars = {x1, x2};//, x3};

   // System's dynamics
   ex dx1 = 0;
   ex dx2 = 0;
   if(!cheb){
   	//minmax approx with deg=3
   	vector< double > coeff {0.9266, -0.1854, 0.0239, -0.001494, 0.00003511};
   	dx1 = dx1 + coeff[0]*pow(x1,1) + coeff[1]*pow(x1,3) + coeff[2]*pow(x1,5) + coeff[3]*pow(x1,7) + coeff[4]*pow(x1,9);
   	dx2 = dx2 + coeff[0]*pow(x2,1) + coeff[1]*pow(x2,3) + coeff[2]*pow(x2,5) + coeff[3]*pow(x2,7) + coeff[4]*pow(x2,9);
   } else {
   	//Chebyshev polynomial with [p,q]=[11,11]
   	vector< double > coeff {0, 1.04549356409021, 0, -0.280020190041206, 0, 0.0613630803596508, 0, -0.00960499014290683, 0, 0.00109138756594997, 0, -9.13173738692450e-05, 0, 5.64415976041858e-06, 0, -2.55223057794851e-07, 0, 8.22486543832501e-09, 0, -1.79184079109511e-10, 0, 2.36842846922030e-12, 0, -1.43600382996810e-14};
   	for(int i=1; i<20; i++){
   		dx1 = dx1 + coeff[i]*pow(x1,i);
		dx2 = dx2 + coeff[i]*pow(x2,i);
   	}
   }
   lst dyns;
   dyns = {dx1,dx2};//,dx3};

   this->vars = vars;
   this->dyns = dyns;


   // Init reach set
   Bundle *B;
   int num_dirs = 4;//3;		// number of bundle directions
   int num_temps = 2;		// number of bundle templates

   // Directions matrix
   vector< double > Li (dim_sys,0);
   vector< vector< double > > L (num_dirs,Li);
   L[0][0] = 1;
   L[1][1] = 1;
   L[2][0] = 1;
   L[2][1] = 1;
   L[3][0] = 1;
   L[3][1] = -1;

   // Template matrix
   vector< int > Ti (dim_sys,0);
   vector< vector< int > > T (num_temps,Ti);
   T[0][0] = 0; T[0][1] = 1;
   T[1][0] = 2; T[1][1] = 3;

   // Offsets for the set of initial conditions
   vector< double > offp (num_dirs,0);
   vector< double > offm (num_dirs,0);
   offp[0] = 1.5; offm[0] = 1.5;
   offp[1] = 1.5; offm[1] = 1.5;
   offp[2] = 2; offm[2] = 2;
   offp[3] = 2; offm[3] = 2;

   B = new Bundle(L,offp,offm,T);

   this->reachSet = B;

}
