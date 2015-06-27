#ifndef SPECKLEAUX_H
#define SPECKLEAUX_H

#include <math.h>
#include <stdio.h>      // fprintf
#include <stdlib.h>     // exit, EXIT_FAILURE

void polint(double* xa, double* ya,
            int n, double x,
            double &y, double &dy);

void polin3(double *x1a,double *x2a,double *x3a,double *ya,
            int l,int m,int n,
            double x1,double x2,double x3,
            double &y, double &dy);

float bessj1(float x);

#endif
