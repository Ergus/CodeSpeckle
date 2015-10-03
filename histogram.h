#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#define dmin(x, y) (((x) < (y)) ? (x) : (y))
#define dmax(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

class speckle;
#include "specklemod.h"

/********************************************//**
 * \brief Usage for the histogram class
 *
 * The histogram class is made to store the values of processing the results
 * The compatibility enforces to make this class as internal for the speckle class.
 * It have to be friend of the speckle (any calculator in general), and should contain
 * a process routine for all the results groups in the speckle class.
 *
 * \author Jimmy Aguilar Mena
 * \version 0.1
 ***********************************************/
class histogram{
    public:
        /// Histogram constructor
        /** \param [in] outter is the outer speckle class. 
            The constructor allocates all the internal arrays and open the
            output file. 
            \param [in] odeltaE bin size for the histogram*/

        histogram(speckle* outter, double odeltaE);

        /// Histogram class destructor
        /** Frees the memory and close the files. */
        ~histogram(){
            free(nos);
            free(countr);            
            free(dos);
            free(sumr);
            free(meanr);
                    
            fclose(f10);
            }

        /// Process routine for colect data in the histograms
        /** \param np [in] Number of elements to process
            \param values [in] Array with the values to process.
            All the data are stored in the arrays during calculations, then 
            the results are written in the output file.
            \note The file is not closed until destructor call. Here only flush 
            and rewind every time. */
        int process(const int np,const double *values);
        
    private:
        const double deltaE;
        int ndeltaE;
        double EMax;
        int nrealiz;
        double *dos, *sumr, *meanr, *diffeigen;
        int *nos, *countr;
        FILE* f10;     // dos.dat
    };

#endif
