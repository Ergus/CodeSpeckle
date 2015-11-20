#include "magma_solver.h"

magma_solver_2stage::magma_solver_2stage(int on, bool ovectors,
                           double omin, double omax,
                           int ongpu):
    solver(on,ovectors,omin,omax),
    uplo(MagmaLower),ngpu(ongpu), 
    iwork(NULL),rwork(NULL),work(NULL){

    // initialize magma system
    magma_init();
    
    jobz = (ovectors?MagmaVec:MagmaNoVec);
    range =(min==max?MagmaRangeAll:MagmaRangeV);
    
    
    //This are local temporal 1 element arrays, we can use
    //a variable and a pointer, but for some reason this is
    //the way it is made in magma examples.
    
    magma_int_t aux_iwork[1];
    double aux_rwork[1];    
    magmaDoubleComplex aux_work[1];

    //First call to get right sizes;

    if(ngpu==1){
        #ifdef DEBUG
        printf("Magma will use 1 gpu\n");
        #endif
        magma_zheevdx_2stage(
                      jobz, range, uplo,
                      n, NULL, n, min, max, 0, 0,
                      &m, NULL,
                      aux_work, -1,
                      aux_rwork, -1,
                      aux_iwork, -1,
                      &info
                      );                          
        }
    else{
        #ifdef DEBUG
        printf("Magma will use %d gpu\n",ngpu);
        #endif        
        magma_zheevdx_2stage_m(
                        ngpu,
                        jobz, range, uplo,
                        n, NULL, n, min, max, 0, 0,
                        &m, NULL,
                        aux_work, -1,
                        aux_rwork, -1,
                        aux_iwork, -1,
                        &info
                        );                                      
        }
    
    if(info!=0){
        fprintf(stderr,"Error magma call for allocation info= %d\n",info);
        exit(EXIT_FAILURE);
        }
    
    //Get sizes
    lwork = (magma_int_t) MAGMA_Z_REAL(aux_work[0]);
    liwork = aux_iwork[0];
    lrwork = (magma_int_t) aux_rwork[0];

    #ifdef DEBUG
    printf("Constructing magma solver_2stage\n");
    #endif // DEBUG        

    }

magma_solver_2stage::~magma_solver_2stage(){

    #ifdef DEBUG
    printf("Destructing magma solver_2stage\n");
    #endif // DEBUG
    
    if(rwork) magma_free_cpu(rwork);
    if(iwork) magma_free_cpu(iwork);
    if(work)  magma_free_cpu(work);
    
    magma_finalize();
    }

int magma_solver_2stage::solve(double complex *oA){

    #ifdef DEBUG
    printf("Solving with magma_2stage solver\n");
    #endif // DEBUG
    
    //Cast for Input Matrix
    magmaDoubleComplex *hA = (magmaDoubleComplex *) oA;

    if(jobz==MagmaVec) oV=oA;  //In magma the vectors are set in the input matrix
    
    // Allocate here, every time the solver is called we allocate and deallocate
    // the arrays, because the other functions will need also some memory.
    // Maybe this will impact performance, but we expect that not too much
    // comparing with the main routine
    // Bellow it frees if allocated, allocate, and check allocation each line
    if( work) free( work); magma_zmalloc_cpu( &work, lwork); dbg_mem( work);
    if(iwork) free(iwork); magma_imalloc_cpu(&iwork,liwork); dbg_mem(iwork);
    if(rwork) free(rwork); magma_dmalloc_cpu(&rwork,lrwork); dbg_mem(rwork);

#ifdef DEBUG
    printf("jobz= %s\n",(jobz==MagmaVec?"MagmaVec":"MagmaNoVec"));
#endif            
    
    if(ngpu==1){
        #ifdef DEBUG
        printf("Calling magma solver with 1 gpu\n");
        #endif        
        magma_zheevdx_2stage(
                      jobz, range, uplo,
                      n, hA, n, min, max, 0, 0,
                      &m, w,
                      work, lwork,
                      rwork, lrwork,
                      iwork, liwork,
                      &info
                      );                          
        }
    else{
        #ifdef DEBUG
        printf("Calling magma solver with %d gpu\n",ngpu);
        #endif
        magma_zheevdx_2stage_m(
                        ngpu,
                        jobz, range, uplo,
                        n, hA, n, min, max, 0, 0,
                        &m, w,
                        work, lwork,
                        rwork, lrwork,
                        iwork, liwork,
                        &info
                        );                                      
        }

    if(info!=0){
        fprintf(stderr,"Error in magma call return: info= %d\n",info);
        return(-1);
        }

    magma_free_cpu(rwork); rwork=NULL;
    magma_free_cpu(iwork); iwork=NULL;
    magma_free_cpu(work);   work=NULL;
    
    #ifdef DEBUG
    printf("Solved with magma_2stage solver\n");
    #endif // DEBUG
    
    return 0;
    }
