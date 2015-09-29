#ifndef MKL_SOLVER_H
#define MKL_SOLVER_H

#include "solver.h"
#include <mkl_types.h>
#include <mkl.h>

/***********************************************//**
 *  \brief Usage for the class mkl_solver.
 *
 *  This solver is always compiled because mkl is a dependency of our code.
 *  This is the slower one, but is included for benchmark purposes.
 *  
 *  \author Jimmy Aguilar Mena
 *  \version 0.1
 ************************************************/
class mkl_solver:public solver{
    public:
        /// Constructor for mkl_solver.
        /** \param [in] n dimension to solve
            \param [in] ovectors bool variable that enables eigen vector calculation.
            \param [in] omin lower limit of interesting values
            \param [in] omax upper limit of interesting values */
        mkl_solver(int n, bool ovectors=false,
                   double omin=0.0, double omax=0.0);

        /// Destructor for mkl_solver        
        ~mkl_solver();

        /// Funtion to calculate eigenvalues and optionaly eigenvectors
        /** To acces the results call the general defined routines for all the solvers
            see: #solver */        
        int solve(double complex* A);
    private:
        int layout;
        char jobz, uplo, range;
        double abstol;
    
    };

#endif //MKL_SOLVER_H
