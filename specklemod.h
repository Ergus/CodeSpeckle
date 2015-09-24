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

#include "solver.h"
#include "mkl_solver.h"

#ifdef UMAGMA
#include "magma_solver.h"
#endif
#ifdef UPLASMA
#include "plasma_solver.h"
#endif

using namespace std;

class speckle {
    public:
        speckle(parser *thepar);
        ~speckle();        
    
    private:
        const int nov,                //Order of interpolation
            Nscatterers,      
            npmax;
        int npxpu,                    //npmax+1, backwar compatibility
            Ntot;                     //Dimension for Matrix A
        
        double *VP,                   //This will have dimmension 3
            *xpos;
        double dxi, dx, min, max, size, vDi, focal, vlambda,
            vDi2, focal2, cofactor, rphase;

        bool nrescale, vectors;       //choice of shift and rescale sum2 speckle
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

        // Solver pointer, will be initialized according to the "solver" value
        class solver* thesolver;
        
    public:
        //Variables for fftspeckle
        double complex Vintensity;
        double complex *A;             //This will have dimension 2
        
        //Here starts the functions, all them are public
        int init(int idseed);
        int ftspeckle();
        int defineA();
        
        //To get the results directly from the solver.
        
        double complex *get_oV(){                     //getter for solution vectors
            return thesolver->get_oV();
            };
        double *get_w(){return thesolver->get_w();};  //getter for solutions values                

        //double EXVT(double xco,double yco,double zco);
        //void correlationspeckle(int idseed);

    };

#endif
