/******************************************************
 * This file is part of the Codespeckle++ distribution Copyright (c) 2015 Jimmy
 * Aguilar Mena.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************/

#ifndef specklemod_h
#define specklemod_h

#include <complex.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "mkl_dfti.h"

#include "base_calculator.h"
#include "histogram.h"

#include "solver.h"
#include "mkl_solver.h"

#ifdef UPLASMA
#include "plasma_solver.h"
#endif

#ifdef UMAGMA
#include "magma_solver.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif


#define frand()((double)rand()/(RAND_MAX)) //random number generator in (0,1)

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

        /// Pointer to the #solver class object.
        /** Possible values are:
            - #plasma_solver
            - #magma_solver
            - #magma_solver_2stage
            - #mkl_solver */
        solver* thesolver;

        /// This is correlation speckle a function to test.
        /** This function should be used during developing time to test the
            potential and other properties.*/
        int correlationspeckle(int idseed);
        
        /// Calculate routine
        /** Makes all the calculations and set the results in the expected 
            direccions to be processed. This is the only function that needs
            to be called intead of #ftspeckle, #init or #defineA. 
            But the others are publoc to provide more flexibility.
            \param [in] idseed Random number generator seed
            \return The return value should be 0, else an error ocurred.*/
        int calculate(int idseed);

        /// Print the values that will be used in this speckle.
        /** \param [in,out] ou to print, default standard output
            \param [in] pre char to put before every line, default '#'*/
        void print(FILE* ou=stdout,char pre='#');

        /// Process function, this needs to be defined is mandatory
        int process_serial(int nvalues, double* array);

        /// This function is a bessel polinomial approximation order 1
        /** This is needed in the correlation routine but can be called 
            from everywhere. It was included here to avoid the creation of
            extra files.*/
        inline double bessj1(double);

        /// This function writes the speckle to a file.
        /** The filename depends of the seed used for the creation of the speckle, 
            to prevent errors.*/
        int writespeckle();

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

        /// Pointer to the init function.
        /** Posiible values are:
            - #init_spherical
            - #init_sum2
            - #init_single
            - #init_shell */
        int (speckle::*pointer_init)(int seed);
        
        histogram* thehist;           ///< Pointer to histogram class, constructed
                                      ///< only if process function is called.        
        
        /// Serial version for the getseed routine
        int getseed_serial(){return( it<=end ? it++ : -1 );}

        /// fprintf specific for f18 file that saves a speckle data for every seed.
        /** This function guarantees that write process is made only if the file 
            is open.
            The FILE objects f18 is null in release mode, but in debug mode
            it is created as dirname/prefix_speckletype_seed.out. */
        void f18_printf(const char * format, ... ){
            #ifdef DEBUG
            if(f18){
                va_list args;
                va_start(args, format);
                vfprintf(f18, format, args);
                va_end(args);
                fflush(f18);
                }
            else{
                fprintf(stderr,"Error: f18 is NULL\n");
                printme();
                }
            #endif
            }

        /** \name Run control variables
            Run system check needs this variables to start and stop
            in this implementation a serial for is used, but this can be reimplemented 
            in the future. */                
        ///\{
        int start,                ///< First seed
            end,                  ///< Last seed
            it;                   ///< This will be like an iterator for getseed
        ///\}
        
        int nov,                  ///< Order of interpolation
            Nscatterers,          ///< Number of scatters
            npmax,
            npxpu,                ///< npmax+1, used for backwar compatibility
            Ntot,                 ///< Dimension for Matrix A
            save_interval;
        
        double *VP,                   ///< This will have dimmension 3
            *xpos;                    ///< array for xpos
        
        double dxi, dx,
            min,                      ///< minimum value for values, initialized to 0
            max,                      ///< maximum value for values, initialized to 0
            size,                     ///< size internal for speckle
            vDi, focal, vlambda,
            vDi2, focal2, cofactor, rphase,
            binsize;                  ///< SIze for the bin in the histogram
        
        /** \name Bool flags
            All are false by default, so should be activated manually.*/
        ///\{
        bool nrescale,                ///< Choice of shift and rescale for sum2
            vectors,                  ///< Option to calculate vector or not
            usegnuplot;               ///< Activates the gnuplot dynamic graphs.
        ///\}
        
    private:

        double complex Vintensity;    ///< Intensity value. Variable for fftspeckle.
        double complex *A,            ///< This will have dimension 2
            *Vk;                      ///< This will have dimension 3

                  ///< Full name for the final file if prefix !null

        string speckletype,           ///< Speckle type name.
            solvername,               ///< Solver name
            fprefix,                  ///< Prefix for output, if "null" no output.
            continuefile;             ///< File to restart calculations, will be readed.

        FILE* f18;                    ///< Output file pointer in debug mode

        int last_seed;
        
        friend class histogram;
        
        void use_option(int opt,const char* thearg);
    };

#endif
