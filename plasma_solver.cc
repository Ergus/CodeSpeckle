#include "plasma_solver.h"

plasma_solver::plasma_solver(int n, bool ovectors,
                           double omin, double omax):
    solver(n,ovectors,omin,omax),
    uplo(PlasmaLower),
    abstol(LAPACKE_dlamch_work('s')){

    // initialize plasma system
    PLASMA_Init(0);
    
    jobz = (ovectors?PlasmaVec:PlasmaNoVec);
    range =(min==max?PlasmaAllVec:PlasmaVec);
    
    if(range==PlasmaAllVec){
        info=PLASMA_Alloc_Workspace_zheevd(n,n,&desc);
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
    }

plasma_solver::~plasma_solver(){

    if(oV) free(oV); oV=NULL;
    
    free(desc);
    PLASMA_Finalize();    
    }

int plasma_solver::solve(double complex *oA){

    //Cast for Input Matrix
    PLASMA_Complex64_t *hA=(PLASMA_Complex64_t *)oA;

    if(range=PlasmaAllVec){
        info=PLASMA_zheevd(jobz, uplo, n,
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

    return 0;
    }
