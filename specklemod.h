#ifndef specklemod_h
#define specklemod_h

#include <complex.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include "mkl_dfti.h"
#include <algorithm>
#include <string.h>
#include <getopt.h>

#define frand()((double)rand()/(RAND_MAX))//random number generator

//This is a debbuger macro for fft functions that needs to
//return 0
#define dbg(x) {                                        \
        if(x!=0){                                       \
            printf("Error: %s returned %d, %s:%d\n",    \
                   #x, x, __FILE__,__LINE__);           \
            exit(EXIT_FAILURE);                         \
        }                                               \
    }

#define dbg_mem(x) {                                                    \
        if(x==NULL){                                                    \
            printf("Error: %s is %d, after allocation in %s:%d\n",      \
                   #x, x, __FILE__,__LINE__);                           \
            exit(EXIT_FAILURE);                                         \
        }                                                               \
    }

//This can be made better but more complication is not needed
extern char Specklename[][10];   //This have to be defined in the source file
enum Speckletype{spherical,sum2,single,shell};
Speckletype atospeckle(const char* type);
const char* speckletostr(Speckletype sp);

using namespace std;

class speckle {
    public:
        speckle(int argc, char** argv);
        ~speckle();        
    
        // variables to define the speckle
        int Nscatterers, npmax;
        double size, vDi, focal, vlambda,
            vDi2, focal2, cofactor, rphase;
        Speckletype speckletype;
        bool nrescale;                 //choice of shift and rescale sum2 speckle
    private:
        int npxpu,
            nov;                       //Order of interpolation
        double *VP,                    //This will have dimmension 3
            *xpos;
        double dxi, dx;
        
        double complex *Vk;

        //Variables for derivative
        double *VPDX, *VPDY, *VPDZ;    //This will have dimmension 3
        double xco, yco, zco;
        
        //private inits for individual cases
        int init_spherical(int idseed);
        int init_sum2(int idseed);
        int init_single(int idseed);
        int init_shell(int idseed);
        
        struct option *longopts;
        int largc; char **largv;       //standard input arguments
        void parser();
        void print_help();
        void use_option(int opt, const char* thearg);

    public:
        //Variables for fftspeckle
        double complex Vintensity;
        //Variables for Plasma/Magma
        int cores, gpus;
        double complex *A;             //This will have dimension 2
        double *w;
        
        //Here starts the functions, all them are public
        int init(int idseed);
        
        double EXVT(double xco,double yco,double zco);
        void correlationspeckle(int idseed);
        int ftspeckle();
        void defineA();

        void print(FILE* output=stdout);
        /*
          void writespeckle();
          void readspeckle();
          void averagecorrelation();
          void diagonalize();
          void DOSspeckle();
        */
    };

#endif
