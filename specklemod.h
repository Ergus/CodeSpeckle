#ifndef specklemod_h
#define specklemod_h

#include <complex.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#define frand()((double)rand()/(RAND_MAX))//random number generator

enum Speckletype{spherical,sum2,single,shell};
extern char Specklename[][10];   //This have to be defined in the source file

class speckle {
    public:
        speckle(int onpmax, Speckletype ospeckletype);
        ~speckle();        
    
        // variables to define the speckle
        const int Nscatterers,
            npmax;
        double size, vDi, focal, vlambda,
               vDi2, focal2, cofactor, rphase;
        const Speckletype speckletype;
        bool nrescale;                 //choice of shift and rescale sum2 speckle        
    private:
        int npxpu;
        double *VP,                    //This will have dimmension 3
               *xpos, *x;
        double dxi, dx;
        const int nov;                 //Order of interpolation
        double complex *Vk;

        //Variables for derivative
        double *VPDX, *VPDY, *VPDZ;    //This will have dimmension 3
        double xco, yco, zco;

    public:
        //Variables for fftspeckle
        double complex Vintensity;
        //Variables for Plasma/Magma
        int cores, gpus;
        double complex *A;             //This will have dimension 2
        double *w;

        //Here starts the functions, all are public
        void init(int idseed);
        /*        double EXVT(double xco,double yco,double zco);
        void correlationspeckle(int idseed);
        void readspeckle();
        void writespeckle();
        void averagecorrelation();
        void ftspeckle();
        void diagonalize();
        void defineA();
        void DOSspeckle();
        */
    };

#endif
