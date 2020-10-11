/**
 * @file Flowpipe.h
 * Represent and manipulate flowpipes of bundle
 * Used to represent the reachable set of a dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef FLOWPIPE_H_
#define FLOWPIPE_H_

#include "Common.h"
#include "Bundle.h"

class Flowpipe {

private:
	vector< Bundle* > flowpipe;			// flowpipe

public:

	// constructors
	Flowpipe();
	Flowpipe(vector< Bundle* >);

	Bundle* get(int i);	// get i-th bundle

	void append( Bundle* bundle );

	int size(){ return this->flowpipe.size(); }

	void print();
	void plotRegion();
	void plotRegionToFile(char *file_name, char color);
	void plotProjToFile(int var, double time_step, char *file_name, char color);

	virtual ~Flowpipe();
};

#endif /* BUNDLE_H_ */
