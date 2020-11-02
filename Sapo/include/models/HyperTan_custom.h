/**
 * @file HyperTan3D.h
 * HYPERTAN activation function
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>, Eleonora.
 * @version 0.1
 */

#ifndef HYPERTAN_CUSTOM_H_
#define HYPERTAN_CUSTOM_H_

#include "Model.h"

class HyperTan_custom : public Model {

private:

public:
	HyperTan_custom(Bundle *B, int dim_sys, double* coeff, int deg);
};

#endif /* HYPERTAN_CUSTOM_H_ */
