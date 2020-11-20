/**
 * @file main.cpp
 * main: This main file reproduces the experiments reported in "Sapo: Reachability Computation and Parameter Synthesis of Polynomial Dynamical Systems"
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include <stdio.h>
#include <iostream>

#include "Common.h"
#include "Bundle.h"
#include "Sapo.h"

#include "VanDerPol.h"
#include "Rossler.h"
#include "SIR.h"
#include "LotkaVolterra.h"
#include "Phosphorelay.h"
#include "Quadcopter.h"
#include "HyperTan3D.h"
#include "HyperTan.h"
#include "HyperTan_custom.h"

#include "SIRp.h"
#include "Influenza.h"
#include "Ebola.h"

using namespace std;

int main(int argc,char** argv){

  // Sapo's options
  sapo_opt options;
  options.trans = 1;			 // Set transformation (0=OFO, 1=AFO)
  options.decomp = 0;			  // Template decomposition (0=no, 1=yes)
  //options.alpha = 0.5;		// Weight for bundle size/orthgonal proximity
  options.verbose = false;

  // Load modles
  vector< Model* > reach_models;
  vector< int > reach_steps;
  // reach_models.push_back(new VanDerPol());      reach_steps.push_back(300);
  // reach_models.push_back(new Rossler());        reach_steps.push_back(250);
  // reach_models.push_back(new SIR(false));       reach_steps.push_back(300);
  // reach_models.push_back(new LotkaVolterra());  reach_steps.push_back(500);
  // reach_models.push_back(new Phosphorelay());   reach_steps.push_back(200);
  // reach_models.push_back(new Quadcopter());     reach_steps.push_back(300);
   int dim_sys = 3;
   int num_dirs = 6;		// number of bundle directions
   int num_temps = 2;		// number of bundle templates

   Bundle *B;
   // Directions matrix
   vector< double > Li (dim_sys,0);
   vector< vector< double > > L (num_dirs,Li);
   L[0][0] = 1;
   L[1][1] = 1;
   L[2][2] = 1;
   L[3][0] = 1;
   L[3][1] = 1;
   L[4][0] = 1;
   L[4][2] = 1;
   L[5][1] = 1;
   L[5][2] = 1;

   // Template matrix
   vector< int > Ti (dim_sys,0);
   vector< vector< int > > T (num_temps,Ti);
   T[0][0] = 0; T[0][1] = 1; T[0][2] = 2;
   T[1][0] = 3; T[1][1] = 4; T[1][2] = 5;

   // Offsets for the set of initial conditions
   vector< double > offp (num_dirs,0);
   vector< double > offm (num_dirs,0);
   offp[0] = 1; offm[0] = 1;
   offp[1] = 1; offm[1] = 0;
   offp[2] = 3; offm[2] = -0.5;
   offp[3] = 1; offm[3] = 0.5;
   offp[4] = 2; offm[4] = -0.5;
   offp[5] = 3; offm[5] = -1;

   B = new Bundle(L,offp,offm,T);

//Nikos
  vector< float > p_coeff {0.9568,0, -0.2107,0,0.02354};
  int deg=p_coeff.size();
  //reach_models.push_back(new HyperTan3D(B, dim_sys));     reach_steps.push_back(1);
  reach_models.push_back(new HyperTan_custom(B, dim_sys, p_coeff, deg));     reach_steps.push_back(1);

  Flowpipe* flowpipe;
  // Compute reach sets
  for(int i=0; i<reach_models.size(); i++){

    cout<<"Model: "<<reach_models[i]->getName()<<"\tReach steps: "<<reach_steps[i]<<"\t";

    Sapo *sapo = new Sapo(reach_models[i],options);
    flowpipe = sapo->reach(reach_models[i]->getReachSet(),reach_steps[i]);	// reachability analysis
  }
  cout<<"\n";

  cout<<"FIGURE Hyperbolic Tangent"<<endl;

  // Generate matlab script to plot flowpipe
  char figTan[] = "plotFigureTan.m";
  flowpipe->plotRegionToFile(figTan,'w');
  // Set picture appearence
  ofstream matlab_script;
  matlab_script.open (figTan, std::ofstream::out | std::ofstream::app);
  matlab_script<<"xlabel('x');\n";
  matlab_script<<"ylabel('y');\n";
  matlab_script<<"axis([-2 2 -2 2 -1 3]);\n";
  // matlab_script<<"view([74 23]);\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<figTan<<" generated\n"<<endl;

  exit(EXIT_SUCCESS);
}
