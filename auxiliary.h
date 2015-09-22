#ifndef AUX_H
#define AUX_H

#include <math.h>
#include <stdio.h>      
#include <stdlib.h>     // exit, EXIT_FAILURE

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

void polint(double* xa, double* ya,
            int n, double x,
            double &y, double &dy);

void polin3(double *x1a,double *x2a,double *x3a,double *ya,
            int l,int m,int n,
            double x1,double x2,double x3,
            double &y, double &dy);

double bessj1(double x);

#endif
