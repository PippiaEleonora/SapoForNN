/**
 * @file HyperTan3D.h
 * HYPERTAN activation function
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>, Eleonora.
 * @version 0.1
 */

#ifndef HYPERTAN3D_H_
#define HYPERTAN3D_H_

#include "Model.h"

class HyperTan3D : public Model {

private:

public:
	HyperTan3D(Bundle *B, int dim_sys);
};

#endif /* HYPERTAN3D_H_ */
