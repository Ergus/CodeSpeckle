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
 * This is an specific class to be used with speckle, so many dependencies are needed.
 * To reuse this clas it is almost needed reimplement it constructor completely.
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

        histogram(speckle* outter);

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

        /// This function reads a prefix_values.dat
        int load_previous(const char *filename){
            STARTDBG
            FILE *fp=fopen(filename,"r");
            if (!fp){
                fprintf(stderr,"Impossible import file %s to continue calculations\n",
                        filename);
                return -1;
                }
            char ignore[1024];
            fgets(ignore, sizeof(ignore), fp);  //ignore first 2 lines
            fgets(ignore, sizeof(ignore), fp);
            //start parsing
            for(int i=0;i<ndeltaE;i++){
                int matched=fscanf(fp,"%*d %*lf %lf %lf %lf %d",
                                   &dos[i], &meanr[i], &sumr[i], &countr[i]);

                // an error exit the program because some
                // values maybe were modified.
                if (matched<4){  
                    fprintf(stderr,"Error importing file %s line %d\n",
                            filename,i+2);
                    exit(EXIT_FAILURE);
                    }
                }
            fclose(fp);
            ENDDBG;
            return(0);
            }
        
    private:
        const double deltaE;
        int ndeltaE;
        double EMax;
        int nrealiz;
        double *dos, *sumr, *meanr, *diffeigen;
        int *nos, *countr;
        FILE* f10;     // dos.dat
        speckle *caller;
    };

#endif
