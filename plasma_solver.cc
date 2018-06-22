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

#include "plasma_solver.h"

plasma_solver::plasma_solver(int n, bool ovectors,
                             double omin, double omax,
                             int ncpu, int start_cpu):
              solver(n,ovectors,omin,omax),
              uplo(PlasmaLower),
              abstol(-1){
    STARTDBG;
    // initialize plasma system
    //omp_set_num_threads(1);
    //goto_set_num_threads(1);
    //mkl_set_num_threads(1);

    coresbind=(int *) malloc(ncpu*sizeof(int));
    for(int i=0;i<ncpu;i++){
        coresbind[i]=start_cpu+i;
        }
    
    PLASMA_Init_Affinity(ncpu,coresbind);
    
    jobz = (ovectors?PlasmaVec:PlasmaNoVec);
    range =(min==max?PlasmaAllVec:PlasmaVec);
    
    if(range==PlasmaAllVec){
        info=PLASMA_Alloc_Workspace_zheev(n,n,&desc);
        }
    else{ //range==PlasmaVec
        info=PLASMA_Alloc_Workspace_zheevr(n,n,&desc);
        }
    
    if(info!=PLASMA_SUCCESS){
        fprintf(stderr,"Error plasma allocation descriptor return %d\n",info);
        printme();
        exit(EXIT_FAILURE);
        }

    if(jobz==PlasmaVec){
        oV=(double complex*) malloc(n*n*sizeof(double complex));
        if(!oV){             //check if the array w was allocated 
            fprintf(stderr,"Error allocating oV (vector for solutions)\n");
            printme();
            exit(EXIT_FAILURE);
            }
        }

    ENDDBG;
    }

plasma_solver::~plasma_solver(){
    STARTDBG;
    
    if(oV) free(oV); oV=NULL;

    free(desc);
    PLASMA_Finalize();
    free(coresbind);
    ENDDBG;
    }

int plasma_solver::solve(double complex *oA){
    STARTDBG;
    
    //Cast for Input Matrix
    PLASMA_Complex64_t *hA=(PLASMA_Complex64_t *)oA;

    if(range==PlasmaAllVec){
        info=PLASMA_zheev(jobz, uplo, n,
                           hA, n, w, desc,
                           oV, n);
        }
    else{
        info=PLASMA_zheevr( jobz, range, uplo,
                            n, hA, n, min, max,
                            0, 0, abstol,
                            &m, w, desc,
                            oV, n );
        }

    if(info!=PLASMA_SUCCESS){
        fprintf(stderr,"Error in Plasma call return: info= %d\n",info);
        return(-1);
        }
    
    ENDDBG;
    return 0;
    }
