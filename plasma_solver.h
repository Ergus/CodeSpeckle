#ifndef PLASMA_SOLVER_H
#define PLASMA_SOLVER_H

#include "solver.h"
#include <plasma.h>
//#include <lapacke.h>

//Here are declared both magma solvers.

class plasma_solver:public solver{
    public:
        plasma_solver(int n, bool ovectors=false,
                      double omin=0.0, double omax=0.0);

        ~plasma_solver();
        
        int solve(double complex* oA);
        
    private:
        PLASMA_enum jobz, uplo, range;
        PLASMA_desc *desc;   //workspace
        int info;            
        const double abstol; // tolerance not used in this version of plasma
    };


#endif //PLASMA_SOLVER_H
