#ifndef SOLVER_H
#define SOLVER_H

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

//==================================================================
// I added the macros here again to allow the usage of the solver in different codes

#ifndef printme
#ifdef MPI_VERSION
#define printme() {                                                  \
        fprintf(stderr,"%s in %s:%d (process %s in %s)\n",           \
                __PRETTY_FUNCTION__, __FILE__, __LINE__,             \
                getenv("OMPI_COMM_WORLD_RANK"), getenv("HOSTNAME")); \
        }
#else  //MPI_VERSION
#define printme() {                                        \
        fprintf(stderr,"%s in %s:%d\n",                    \
                __PRETTY_FUNCTION__, __FILE__, __LINE__);  \
        }
#endif //MPI_VERSION
#endif // printme

#ifdef DEBUG;
 #define STARTDBG fprintf(stderr,"Calling %s: %s:%d\n",                    \
                         __PRETTY_FUNCTION__,__FILE__,__LINE__);
 #define ENDDBG fprintf(stderr,"Ended %s: %s:%d\n",                        \
                         __PRETTY_FUNCTION__,__FILE__,__LINE__);
#else
 #define STARTDBG
 #define ENDDBG
#endif

#ifndef dbg
#define dbg(x) {                                                           \
        if(x!=0){                                                          \
            fprintf(stderr,"Error: %s returned %ld\n",#x, x);              \
            printme();                                                     \
            return(-1);                                                    \
            }                                                              \
        }
#endif  //dbg

#ifndef dbg_mem
#define dbg_mem(x) {                                                       \
        if(x==NULL){                                                       \
            fprintf(stderr,"Error: %s is %p, after allocation\n",#x, x);   \
            printme();                                                     \
            return(-1);                                                    \
            }                                                              \
        }
#endif  //dbg_mem
//==================================================================

/***********************************************//**
 *  \brief Usage for the virtual class solver.
 *
 *  This is the virtual class "solver", all the possible solvers implemented and 
 *  that will be implemented in the future should therive from this base class.
 *  This represents an abstract class that works as interface betwen the speckle 
 *  and the solver libraries.
 *  The solvers should be constructed for specific sizes.
 *  As there are different interfaces the interface works in this way:
 *  constructs ones, pass all the needed parameters, then call solve(A)
 *  The number of values are in m (int get_m()), the values in w (double* get_w)
 *  and the vectors if calculated in (double* get_oV) 
 *  
 *  \note Modifications here will affect all the solvers.
 *  \author Jimmy Aguilar Mena
 *  \version 0.1
 ************************************************/
class solver{
    public:
        /// Virtual solve method
        /** The virtual intance for a general solve methos that should be implemented
            in every different derived class
            \param [in] oA Matrix to diagolalize. The size should be fixed in the 
            constructor of the solver */
        virtual int solve(double complex *oA) = 0; 

        /// Getter for number of solutions
        int get_m(){return m;};

        /// Getter for solutions values
        double *get_w(){return w;};
        
        /// Getter vectors (if calculated)
        double complex *get_oV(){return oV;};

        /// Destructor for the solver.
        virtual ~solver(){
            if(w) free(w); w=NULL;
            }        
        
    protected:
        /// Virtual constructor for the solver
        /** \param on Matrix dimension to solve
            \param [in] ovectors bool variable to indicate if solve eigenvectors.
            - true: Solve also eigenvectors
            - false: Solve only eigenvalues.
            The interface internally should care about the use of the right routine
            \param [in] omin lower bound of the interesting interval
            \param [in] omax upper bound of the interesting interval
            \note: if min==max then all the spectrum is computed */
        solver(int on, bool ovectors, double omin, double omax):
            n(on),              // Matrix Dimension
            vectors(ovectors),  // Calculate also vectors or not
            min(omin),          // min value in range
            max(omax),          // max value in range
            m(on),              // initialize m=n to the max possible value 
            oV(NULL){           // this returns the number of values found
            w=(double*) malloc(on*sizeof(double));
            if(!w){             //check if the array w was allocated 
                fprintf(stderr,"Error allocating w (vector for solutions)\n");
                printme();
                exit(EXIT_FAILURE);
                }
            }
        
        const int n;         ///< Dimension for the matrix to solve (internal storage).
        const bool vectors;  ///< Remember if it should calculate vectors or not
        /// Internalstorage for interval bounds
        ///\{
        double min, max;
        ///\}
        int m,               ///< This will remember amount of values found
                             ///< in last run inside the interesting interval
            info;            ///< Used for error checking, gets the output of                                        ///< library calls to check errors
        
        double *w;           ///< Array for eigen values
        double complex *oV;  ///< Array for eigen vectors
    };


#endif //SOLVER_H
