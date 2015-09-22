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
#include "auxiliary.h"

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
