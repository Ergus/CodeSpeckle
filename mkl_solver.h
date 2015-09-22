#ifndef MKL_SOLVER_H
#define MKL_SOLVER_H

#include "solver.h"
#include <mkl_types.h>
#include <mkl.h>

// Declaration of MKL solver class.
// in parent class solver are deffined n,
// vectors, min, max, m, w, oV

class mkl_solver:public solver{
    public:
        mkl_solver(int n, bool ovectors=false,
                   double omin=0.0, double omax=0.0);

        ~mkl_solver();
        
        int solve(double complex* A);
    private:
        int layout;
        char jobz, uplo, range;
        double abstol;
    
    };

#endif //MKL_SOLVER_H
