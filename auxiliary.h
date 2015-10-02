#ifndef AUX_H
#define AUX_H

#include <math.h>
#include <stdio.h>      
#include <stdlib.h>     // exit, EXIT_FAILURE
#include <stdarg.h>     // argument list dbg_printf

#define frand()((double)rand()/(RAND_MAX)) //random number generator in (0,1)

//This is a debbuger macro for fft functions that needs to
//return 0 each speckle will print its own error message before interrupt

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
 #define STARTDBG fprintf(stderr,"Calling %s: %s:%d\n",                   \
                         __PRETTY_FUNCTION__,__FILE__,__LINE__);
 #define ENDDBG fprintf(stderr,"Ended %s: %s:%d\n",                       \
                         __PRETTY_FUNCTION__,__FILE__,__LINE__);
#else
 #define STARTDBG
 #define ENDDBG
#endif

#ifndef dbg
#define dbg(x) {                                                           \
        if(x!=0){                                                          \
            fprintf(stderr,"Error: %s returned %d\n",#x, x);               \
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

/// A printf implementation that only prints in debug mode.
/** The arguments are the same than printf, but no information 
will be printed in starndard outpur if the DEBUG macro is not 
defined at compile time*/
void dbg_printf(const char * format, ... );

/// Polinomial interpolation 1 dimension.
void polint(double* xa, double* ya,
            int n, double x,
            double &y, double &dy);

/// Polinomial interpolation in 3 dimensions
void polin3(double *x1a,double *x2a,double *x3a,double *ya,
            int l,int m,int n,
            double x1,double x2,double x3,
            double &y, double &dy);

double bessj1(double x);

#endif
