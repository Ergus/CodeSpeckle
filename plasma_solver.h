#ifndef PLASMA_SOLVER_H
#define PLASMA_SOLVER_H

#include "solver.h"
#include <plasma.h>

/***********************************************//**
 *  \brief Usage for the class plasma_solver.
 *
 *  This solver is compiled only if the PLASMA enviroment variables are set.
 *  Else all the references to it are eliminated in compile time.
 *  As all the solvers in the code the base class for this in the solver class.
 *  
 *  \author Jimmy Aguilar Mena
 *  \version 0.1
 ************************************************/
class plasma_solver:public solver{
    public:
        /// Constructor for plasma_solver.
        /** \param [in] n dimension to solve
            \param [in] ovectors bool variable that enables eigen vector calculation.
            \param [in] omin lower limit of interesting values
            \param [in] omax upper limit of interesting values
            \param [in] ncpu number of cored to be used in plasma execution
            \param [in] start_cpu index for the start core for the process, the default value is zero, but mpi interface needs to set it for a good affinity.*/
        plasma_solver(int n, bool ovectors=false,
                      double omin=0.0, double omax=0.0,
                      int ncpu=0, int start_cpu=0);
        
        /// Destructor for plasma_solver
        ~plasma_solver();

        /// Funtion to calculate eigenvalues and optionaly eigenvectors
        /** To acces the results call the general defined routines for all the solvers
            see: #solver */
        int solve(double complex* oA);
        
    private:
        PLASMA_enum jobz, uplo, range;
        PLASMA_desc *desc;   //workspace
        int info;            
        const double abstol; // tolerance not used in this version of plasma
        int *coresbind;
    };


#endif //PLASMA_SOLVER_H
