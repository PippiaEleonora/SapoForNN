#include <stdio.h>
#include <iostream>

#include "Common.h"
#include "Bundle.h"
#include "Sapo.h"

#include "HyperTan3D.h"


using namespace std;


extern "C" {
	int computeSapo(int dim_sys, int num_dirs, int num_temps, double **cL, double **cT, double *coffp, double *coffm, double **A){
		// Sapo's options

  sapo_opt options;
  options.trans = 1;			 // Set transformation (0=OFO, 1=AFO)
  options.decomp = 0;			  // Template decomposition (0=no, 1=yes)
  //options.alpha = 0.5;		// Weight for bundle size/orthgonal proximity
  options.verbose = false;

  // Read input values
  // 1) direction matrix
  vector< double > Li (dim_sys,0);
  vector< vector< double > > L (num_dirs,Li);
  for(int i=0; i<num_dirs; i++){
	for(int j=0; j<dim_sys; j++){
		L[i][j] = cL[i][j];
	}
  }

  // 2) template matrix
  vector< int > Ti (dim_sys,0);
  vector< vector< int > > T (num_temps,Ti);
  for(int i=0; i<num_temps; i++){
	for(int j=0; j<dim_sys; j++){
		T[i][j] = cT[i][j];
	}
  }

  // 3) offsets for the set of initial conditions
  vector< double > offp (num_dirs,0);
  vector< double > offm (num_dirs,0);
  for(int i=0; i<num_dirs; i++){
	offp[i] = coffp[i];
	offm[i] = coffm[i];
  }

  Bundle *B;
  B = new Bundle(L,offp,offm,T);

  // Load modles
  Model* reach_model = new HyperTan3D(B, dim_sys);
  int reach_steps = 1;

  // Compute reach sets
  Flowpipe* flowpipe;
  cout<<"Model: "<<reach_model->getName()<<"\tReach steps: "<<reach_steps<<"\t";
  Sapo *sapo = new Sapo(reach_model,options);
  flowpipe = sapo->reach(reach_model->getReachSet(),reach_steps);	// reachability analysis
	
  // Create the output matrix
  Bundle *Bout;
  Bout = flowpipe->get(1);

  vector<double> offp_out = Bout->getoffps();
  vector<double> offm_out = Bout->getoffms();
  vector<vector<double>> L_out = Bout->getLs();
  int n_cons = Bout->getSize();

  for(int i=0; i<n_cons; i++){
	A[i][0] = offp_out[i];
	for(int j=0; j<dim_sys; j++){
		A[i][j+1] = L_out[i][j];
	}
  }
  for(int i=0; i<n_cons; i++){
	A[n_cons+i][0] = offm_out[i];
	for(int j=0; j<dim_sys; j++){
		A[n_cons+i][j+1] = -L_out[i][j];
	}
  }

  /*cout<<bp[0]<<endl;
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
  cout<<figTan<<" generated\n"<<endl;*/


  return 2*n_cons;
	}
}
