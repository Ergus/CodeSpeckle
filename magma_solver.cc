#include "magma_solver.h"

magma_solver::magma_solver(int on, bool ovectors,
                           double omin, double omax,
                           int ongpu):
    solver(on,ovectors,omin,omax),
    ngpu(ongpu), uplo(MagmaLower),
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
    if(range==MagmaRangeAll){
        if(ngpu==1){
            magma_zheevd(
                         jobz, uplo, n, NULL, n, NULL,
                         aux_work, -1,
                         aux_rwork, -1,
                         aux_iwork, -1,
                         &info
                         );
            }
        else{
            magma_zheevd_m(
                           ngpu,
                           jobz, uplo, n, NULL, n, NULL,
                           aux_work, -1,
                           aux_rwork, -1,
                           aux_iwork, -1,
                           &info
                           );
            }
        }
    else{
        if(ngpu==1){
            magma_zheevdx(
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
            magma_zheevdx_m(
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
        }
    
    if(info!=0){
        fprintf(stderr,"Error magma call for allocation info= %d\n",info);
        exit(EXIT_FAILURE);
        }
    
    //Get sizes
    lwork = (magma_int_t) MAGMA_Z_REAL(aux_work[0]);
    liwork = aux_iwork[0];
    lrwork = (magma_int_t) aux_rwork[0];

    }

magma_solver::~magma_solver(){
    if(rwork) magma_free_cpu(rwork);
    if(iwork) magma_free_cpu(iwork);
    if(work)  magma_free_cpu(work);
    
    magma_finalize();
    }

int magma_solver::solve(double complex *oA){
    
    //Cast for Input Matrix
    magmaDoubleComplex *hA = (magmaDoubleComplex *) oA;

    if(jobz==MagmaVec) oV=oA;  //In magma the vectors are set in the input matrix
    
    // Allocate here, every time the solver is called we allocate and deallocate
    // the arrays, because the other functions will need also some memory.
    // Maybe this will impact performance, but we expect that not too much
    // comparing with the main routine
    // Bellow frees if allocated, allocate, and check allocation each line
    if( work) free( work); magma_zmalloc_cpu( &work, lwork); dbg_mem( work);
    if(iwork) free(iwork); magma_imalloc_cpu(&iwork,liwork); dbg_mem(iwork);
    if(rwork) free(rwork); magma_dmalloc_cpu(&rwork,lrwork); dbg_mem(rwork);

    if(range==MagmaRangeAll){
        if(ngpu==1){
            magma_zheevd(
                         jobz, uplo,
                         n, hA, n,
                         w,
                         work, lwork,
                         rwork, lrwork,
                         iwork, liwork,
                         &info
                         );
            }
        else{
            magma_zheevd_m(
                           ngpu,
                           jobz, uplo,
                           n, hA, n,
                           w,
                           work, lwork,
                           rwork, lrwork,
                           iwork, liwork,
                           &info
                           );
            }

        }
    else{
        if(ngpu==1){
            magma_zheevdx(
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
            magma_zheevdx_m(
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
        }

    if(info!=0){
        fprintf(stderr,"Error in magma call return: info= %d\n",info);
        return(-1);
        }

    magma_free_cpu(rwork); rwork=NULL;
    magma_free_cpu(iwork); iwork=NULL;
    magma_free_cpu(work);   work=NULL;

    return 0;
    }
