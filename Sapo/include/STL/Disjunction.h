/**
 * @file Disjunction.h
 * Disjunction STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef DISJUNCTION_H_
#define DISJUNCTION_H_

#include "STL.h"

class Disjunction : public STL {

private:
	STL * f1, * f2;	// subformulas

public:

	Disjunction(STL * f1, STL * f2);

	STL * getLeftSubFormula(){return f1;};
	STL * getRightSubFormula(){return f2;};

	void print();

	virtual ~Disjunction();
};

#endif /* CONJUNCTION_H_ */
