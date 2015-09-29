#ifndef specklemod_h
#define specklemod_h

#include <complex.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include "mkl_dfti.h"
#include <string.h>
#include <string>
#include <getopt.h>

#include "auxiliary.h"

#include "base_calculator.h"
#include "histogram.h"

#include "solver.h"
#include "mkl_solver.h"

#define DEFAULT_SOLVER mkl_solver

#ifdef UPLASMA
#include "plasma_solver.h"
#undef DEFAULT_SOLVER
#define DEFAULT_SOLVER plasma_solver
#endif

#ifdef UMAGMA
#include "magma_solver.h"
#undef DEFAULT_SOLVER
#define DEFAULT_SOLVER magma_solver  //the default is not 2stage version
#endif

using namespace std;
class histogram;
/********************************************//**
 * \brief Usage for the speckle class.
 *
 * This is the speckle class. It contains the init routines, fftw and calls to
 * diagonalize routines.
 * Values are stored to be sended and procesed in the master node.
 * In this moments the solver can calculate also eigenvectors, but no support
 * have been added and no processing function have been defined, so no used for 
 * the moment.
 * 
 * \author Jimmy Aguilar Mena
 * \version 0.1
 ***********************************************/
class speckle: public  base_calculator{
    public:
        /// Speckle basic constructor.
        /** \param [in] argc standard input counter
            \param [in] argv standard input char array
            The function will call a parser routine for the input file and some 
            command line arguments.*/
        speckle(int argc, char **argv);

        /// Speckle destructor.
        /** It frees the memory arrays, deletes the solver and clode the files.*/
        ~speckle();
        
        /// Speckle initialization sub-routine
        /** This will call the right init routine throw the pointer 
            pointer_init initialized in the constructor.
            \param idseed Random number generator seed
            \return The return value should be 0, else an error ocurred.*/
        int init(int idseed);

        /// Fast Fourier transform routine
        /** This routine calls the fft implementations in mkl. 
            And generates the transform in Vk.*/
        int ftspeckle();
        
        /// Definition for the A array to be diagonalized.
        int defineA();

        /// Calculate routine
        /** Makes all the calculations and set the results in the expected 
            direccions to be processed
            \param [in] idseed Random number generator seed
            \return The return value should be 0, else an error ocurred.*/
        int calculate(int idseed){
            STARTDBG
            dbg(init(idseed));
            dbg(ftspeckle());
            dbg(defineA());
            dbg(thesolver->solve(A));
            indices=thesolver->get_m();
            values=thesolver->get_w();
            ENDDBG
            return 0;
            }        

        /// Print the values that will be used in this speckle.
        /** \param [in,out] ou to print, default standard output
            \param [in] pre char to put before every line, default '#'*/
        void print(FILE* ou=stdout,char pre='#');

        /// Process function, this needs to be defined is mandatory
        virtual int process(int nvalues, double* array);

    protected:
        /** \name Init routines
            This routines initializes the speckle. The #pointer_init points to 
            one of this and is set in the constructor
            \param [in] idseed seed for random number generator
            \return The return value should be 0, else an error ocurred.*/
        ///\{
        int init_spherical(int idseed);    ///< Init for spherica speckle
        int init_sum2(int idseed);         ///< Init for sum2 speckle
        int init_single(int idseed);       ///< Init for single speckle
        int init_shell(int idseed);        ///< Init for shell speckle
        ///\}

        /// Parses input file and command line
        /** The data to parse is suposed to be in largc, and largv. 
            This is called in the constructor.
            \return The return value should be 0, else an error ocurred.*/
        int parse();

        /// Pointer to the #solver class object.
        /** Possible values are:
            - #plasma_solver
            - #magma_solver
            - #magma_solver_2stage
            - #mkl_solver */
        solver* thesolver;

        /// Pointer to the init function.
        /** Posiible values are:
            - #init_spherical
            - #init_sum2
            - #init_single
            - #init_shell */
        int (speckle::*pointer_init)(int seed);

        /// Serial version for the getseed routine
        virtual int getseed(){return( it<=end ? it++ : -1 );}

        histogram* thehist;         ///< Pointer to histogram class, constructed
                                      ///< only if process function is called.        

    protected:
        int nov,                ///< Order of interpolation
            Nscatterers,              ///< Number of scatters
            npmax;
        int npxpu,                    ///< npmax+1, used for backwar compatibility
            Ntot,                     ///< Dimension for Matrix A
            it;
        
        double *VP,                   ///< This will have dimmension 3
            *xpos;                    ///< array for xpos
        
        double dxi, dx,
            micron,
            min,                      ///< minimum value for values, initialized to 0
            max,                      ///< maximum value for values, initialized to 0
            size,                     ///< size internal for speckle
            vDi, focal, vlambda,
            vDi2, focal2, cofactor, rphase;

        bool nrescale,                ///< Choice of shift and rescale for sum2
            vectors;                  ///< Option to calculate vector or not
        double complex *Vk;

        int largc;                    ///< Standard input counter
        char **largv;                 ///< Standard input array

        char outputname[50];          ///< Full name for the final file if not null
        FILE* f18;                    ///< Output file pointer.

        string speckletype,           ///< Speckle type name.
            solver,                   ///< Solver name
            fprefix;                  ///< Prefix for output, if "null" no output.

        double complex Vintensity;    ///< Intensity value. Variable for fftspeckle.
        double complex *A;            ///< This will have dimension 2, matrix for solver
        friend class histogram;
        void use_option(int opt,const char* thearg);

    };

#endif
