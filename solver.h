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

#ifndef dbg
#define dbg(x) {                                                           \
        if(x!=0){                                                          \
            fprintf(stderr,"Error: %s returned %ld\n",#x, x);                 \
            printme();                                                     \
            return(-1);                                                    \
            }                                                              \
        }
#endif  //dbg

#ifndef dbg_mem
#define dbg_mem(x) {                                                       \
        if(x==NULL){                                                       \
            fprintf(stderr,"Error: %s is %p, after allocation\n",#x, x);     \
            printme();                                                     \
            return(-1);                                                    \
            }                                                              \
        }
#endif  //dbg_mem
//==================================================================

//...oooOOO000OOOooo......Solver virtual class......oooOOO000OOOooo...
//...oooOOO0 Modifications here will affect all the solvers 0OOOooo...

// As there are different interfaces the interface works in this way:
// constructs ones, pass all the needed parameters, then call solve(A)
// The number of values are in m (int get_m()), the values in w (double* get_w)
// and the vectors if calculated in (double* get_oV)

class solver{
    public:
        virtual int solve(double complex *oA) = 0; //Hermitic Matrix
        
        int get_m(){return m;};                    //getter for number of solutions
        double *get_w(){return w;};                //getter for solutions values
        double complex *get_oV(){return oV;};      //getter for solution vectors
        virtual ~solver(){
            if(w) free(w); w=NULL;
            }        
        
    protected:
        solver(int on, bool ovectors, double omin, double omax):
            n(on),            // Matrix Dimension
            vectors(ovectors),  // Calculate also vectors or not
            min(omin),          // min value in range
            max(omax),          // max value in range
            m(on),            // initialize m=n to the max possible value 
            oV(NULL){           // this returns the number of values found
            w=(double*) malloc(on*sizeof(double));
            if(!w){             //check if the array w was allocated 
                fprintf(stderr,"Error allocating w (vector for solutions)\n");
                printme();
                exit(EXIT_FAILURE);
                }
            }
        
        const int n;         // Dimension for the matrix.
        const bool vectors;  // Calculate also vectors or not??
        double min, max;
        int m, info;         // This variable will remember the values found in last run

        double *w;           // Vector for eigen values
        double complex *oV;  // Array for eigen vectors
    };


#endif //SOLVER_H
