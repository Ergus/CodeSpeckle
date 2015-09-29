#ifndef MAGMA_SOLVER_H
#define MAGMA_SOLVER_H

#include "solver.h"
#include "magma.h"

/***********************************************//**
 *  \brief Usage for the class magma_solver.
 *
 *  This solver is compiled only if the Magma enviroment variables are set.
 *  Else all the references to it are eliminated in compile time.
 *  As all the solvers in the code the base class for this in the solver class.
 *  This solver is slower than the 2stage one, but requires less memory, that's why 
 *  it is defined as the default.
 *
 *  \author Jimmy Aguilar Mena
 *  \version 0.1
 ************************************************/
class magma_solver:public solver{
    public:
        /// Constructor for magma_solver.
        /** \param [in] on dimension to solve
            \param [in] ovectors bool variable that enables eigen vector calculation.
            \param [in] omin lower limit of interesting values
            \param [in] omax upper limit of interesting values 
            \param [in] ongpu number of gpus to use. */
        magma_solver(int on, bool ovectors=false,
                     double omin=0.0, double omax=0.0,
                     int ongpu=1);
        
        /// Destructor for magma_solver
        ~magma_solver();

        /// Funtion to calculate eigenvalues and optionaly eigenvectors
        /** To acces the results call the general defined routines for all the solvers
            see: #solver */        
        int solve(double complex* oA);
        
    private:
        magma_vec_t jobz;
        const magma_uplo_t uplo;
        magma_range_t range;    
        const magma_int_t ngpu;
        magma_int_t lwork, liwork, lrwork, info;

        //Needed arrays (work space)
        magma_int_t *iwork;
        double *rwork;
        magmaDoubleComplex *work;
    };

/***********************************************//**
 *  \brief Usage for the class magma_solver_2stage.
 *
 *  This solver is compiled only if the Magma enviroment variables are set.
 *  Else all the references to it are eliminated in compile time.
 *  As all the solvers in the code the base class for this in the solver class.
 *  This solver is faster, but requires more memory than de default magma_solver, 
 *  that's why you should specify it on run time.
 *
 *  \author Jimmy Aguilar Mena
 *  \version 0.1
 ************************************************/
class magma_solver_2stage:public solver{
    public:
        /// Constructor for magma_solver_2stage.
        /** \param [in] on dimension to solve
            \param [in] ovectors bool variable that enables eigen vector calculation.
            \param [in] omin lower limit of interesting values
            \param [in] omax upper limit of interesting values 
            \param [in] ongpu number of gpus to use. */
        magma_solver_2stage(int on, bool ovectors=false,
                            double omin=0.0, double omax=0.0,
                            int ongpu=1);

        /// Destructor for magma_solver_2stage
        ~magma_solver_2stage();

        /// Funtion to calculate eigenvalues and optionaly eigenvectors
        /** To acces the results call the general defined routines for all the solvers
            see: #solver */                
        int solve(double complex* oA);
        
    private:
        magma_vec_t jobz;
        const magma_uplo_t uplo;
        magma_range_t range;    
        const magma_int_t ngpu;
        magma_int_t lwork, liwork, lrwork, info;

        //Needed arrays (work space)
        magma_int_t *iwork;
        double *rwork;    
        magmaDoubleComplex *work;    
    };


#endif //MAGMA_SOLVER_H
