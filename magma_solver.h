#ifndef MAGMA_SOLVER_H
#define MAGMA_SOLVER_H

#include "solver.h"
#include "magma.h"
//#include "magma_lapack.h"

//Here are declared both magma solvers.

class magma_solver:public solver{
    public:
        magma_solver(int on, bool ovectors=false,
                     double omin=0.0, double omax=0.0,
                     int ongpu=1);
        ~magma_solver();
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

class magma_solver_2stage:public solver{
    public:
        magma_solver_2stage(int on, bool ovectors=false,
                            double omin=0.0, double omax=0.0,
                            int ongpu=1);
        ~magma_solver_2stage();
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
