/******************************************************
 * This file is part of the Codespeckle++ distribution Copyright (c) 2015 Jimmy
 * Aguilar Mena.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************/

#include "mkl_solver.h"

mkl_solver::mkl_solver(int n, bool ovectors,
                       double omin, double omax, int ncpu):
    solver(n,ovectors,omin,omax),
    layout(LAPACK_ROW_MAJOR),
    uplo('L'){
    STARTDBG
    jobz = (ovectors?'V':'N');   // Values and or vectors
    range =(min==max?'A':'V');   // all or a range
    abstol = -1.0;               // error, <0 means default value
    

    // if it will calculate a range of values and vectors the zheevr 
    // writes values in a different array not as zheevd
    if((jobz=='V')&&(range=='V')){

        mkl_set_num_threads(ncpu);
        
        oV=(double complex *) malloc(n*n*sizeof(double complex));
        if (!oV){
            fprintf(stderr,"Error allocating oV (array for vectors)\n");
            printme();
            exit(EXIT_FAILURE);
            }
        }
    ENDDBG
    }

mkl_solver::~mkl_solver(){
    STARTDBG
    if((jobz=='V')&&(range=='V')&&(oV)) free(oV); oV=NULL;
    
    #ifdef DEBUG
    printf("Destructing mkl solver\n");
    #endif // DEBUG
    ENDDBG
    }

int mkl_solver::solve(double complex *oA){
    STARTDBG    
    //Cast for Input Matrix and output
    MKL_Complex16 *hA=(MKL_Complex16*) oA;
    MKL_Complex16 *toV=(MKL_Complex16*) oV; 
    if(range=='A'){
        // oV points to oA for this function, the other one is different
        if(jobz=='V') oV=oA;
        
        info=LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, uplo, n, hA, n, w);
        if(info!=0){
            fprintf(stderr,"Error in mkl zheev call return: info= %d\n",info);
            printme();
            return(-1);
            }                
        }
    else{ //if range=='V'
        info = LAPACKE_zheevr( layout, jobz, range, uplo, n, hA, n,
                               min, max, 0, 0, abstol, &m, w, toV, n, NULL );
        if(info!=0){
            fprintf(stderr,"Error in mkl zheevr call return: info= %d\n",info);
            printme();
            return(-1);
            }        
        }

    ENDDBG
    return 0;
    }
