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

#include "parser.h"
#include "speckleauxiliary.h"

#define frand()((double)rand()/(RAND_MAX)) //random number generator in (0,1)

//This is a debbuger macro for fft functions that needs to
//return 0 each speckle will print its own error message before interrupt
#ifdef MPI_VERSION
#define printme() {                                                  \
        fprintf(stderr,"%s in %s:%d (process %s in %s)\n",           \
                __PRETTY_FUNCTION__, __FILE__, __LINE__,             \
                getenv("OMPI_COMM_WORLD_RANK"), getenv("HOSTNAME")); \
        }
#else
#define printme() {                                        \
        fprintf(stderr,"%s in %s:%d\n",                    \
                __PRETTY_FUNCTION__, __FILE__, __LINE__);  \
        }
#endif


#define dbg(x) {                                                           \
        if(x!=0){                                                          \
            fprintf(stderr,"Error: %s returned %ld",#x, x);                 \
            printme();                                                     \
            return(-1);                                                    \
            }                                                              \
        }

#define dbg_mem(x) {                                                       \
        if(x==NULL){                                                       \
            fprintf(stderr,"Error: %s is %p, after allocation",#x, x);     \
            printme();                                                     \
            return(-1);                                                    \
            }                                                              \
        }

using namespace std;

class parser;

class speckle {
    public:
        speckle(parser *thepar);
        ~speckle();        
    
        // variables to define the speckle
        int Nscatterers, npmax;
        double size, vDi, focal, vlambda,
            vDi2, focal2, cofactor, rphase;
        bool nrescale;                 //choice of shift and rescale sum2 speckle
    private:
        int npxpu,
            nov;                       //Order of interpolation
        double *VP,                    //This will have dimmension 3
            *xpos;
        double dxi, dx;
        
        double complex *Vk;

        //private inits for individual cases
        int init_spherical(int idseed);
        int init_sum2(int idseed);
        int init_single(int idseed);
        int init_shell(int idseed);
        
        struct option *longopts;
        int largc; char **largv;       //standard input arguments

        int (speckle::*pointer_init)(int seed);
        parser *parptr;                //Parser pointer
        char outputname[50];
        FILE* f18;
        string speckletype, solver, fprefix;
    public:
        //Variables for fftspeckle
        double complex Vintensity;
        double complex *A;             //This will have dimension 2
        double *w;
        
        //Here starts the functions, all them are public
        int init(int idseed);
        int ftspeckle();
        int defineA();

        //double EXVT(double xco,double yco,double zco);
        //void correlationspeckle(int idseed);

    };

#endif
