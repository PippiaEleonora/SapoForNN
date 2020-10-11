/**
 * @file Conjunction.cpp
 * Conjunction STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Conjunction.h"


/**
 * Constructor that instantiates a Conjunction STL formula (f1 /\ f2)
 *
 * @param[in] f1 left conjunct
 * @param[in] f2 right conjunct
 */
Conjunction::Conjunction(STL * f1, STL * f2){
	this->f1=f1;
	this->f2=f2;
	type=CONJUNCTION;
};

/**
 * Print the formula
 */
void Conjunction::print(){
	cout<<"(";
	this->f1->print();
	cout<<") and (";
	this->f2->print();
	cout<<")";
}

Conjunction::~Conjunction() {}
