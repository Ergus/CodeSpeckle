#include "histogram.h"

histogram::histogram(speckle *outter, double odeltaE):
    deltaE(odeltaE),nrealiz(0),
    nos(NULL),countr(NULL),dos(NULL),sumr(NULL),meanr(NULL)
{
    STARTDBG
    const double tmp =2.0*M_PI/outter->size;
    const double tmp2=outter->npmax/2.0;
    EMax=3.*tmp*tmp*tmp2*tmp2+5.0*creal(outter->Vintensity);
		
    ndeltaE=EMax/deltaE+1;
    //int arrays
    nos=(int*) calloc(ndeltaE,sizeof(int));
    countr=(int*) calloc(ndeltaE,sizeof(int));
    //double arrays
    dos=(double*) calloc(ndeltaE,sizeof(double));
    sumr=(double*) calloc(ndeltaE,sizeof(double));
    meanr=(double*) calloc(ndeltaE,sizeof(double));

    if (!(nos && countr && dos && sumr && meanr)){
            fprintf(stderr,"Error: allocating memory for histogram\n");
            printme();
            exit(EXIT_FAILURE);
            }
    
    f10=fopen((outter->fprefix+"_dos.dat").c_str(),"w");       // dos.dat
    ENDDBG
    }

int histogram::process(const int np,const double *values){
    STARTDBG
    nrealiz++;

    double diffeigen, diffeigen_p;
    // Check right border
    int neigen=values[np-1]/deltaE;
    if(neigen<ndeltaE) nos[neigen]++;
    // check first element and decides if starts
    neigen=values[0]/deltaE;
    if(neigen<ndeltaE){
        nos[neigen]++;
        diffeigen_p=values[1]-values[0];
        for(int i=1;i<np-1;i++){
            neigen=values[i]/deltaE;
            if(neigen>=ndeltaE) break;
            nos[neigen]++;            
            diffeigen=diffeigen_p;
            diffeigen_p=values[i+1]-values[i];
            sumr[neigen]+=(dmin(diffeigen,diffeigen_p)
                           /dmax(diffeigen,diffeigen_p));
            countr[neigen]++;
            }
        }
    else{
        fprintf(stderr,"Warning no values in the interesting range\n");
        }
    for(int i=0;i<ndeltaE;i++) dos[i]=((double)nos[i])/nrealiz;

    rewind(f10);
    fprintf(f10,"#Emax= %lf, deltaE= %lf, nrealiz=%d\n",
            EMax,deltaE,nrealiz);
    fprintf(f10,"#n E(n) dos(E) meanr(E) sumr(E) countr(E)\n");
                    
    double E=0;
    for(int i=0;i<ndeltaE;i++){
        E+=deltaE;
        meanr[i]=(countr[i]==0 ? 0.0 : sumr[i]/countr[i]);
        fprintf(f10,"%d %lf %lf %lf %lf %d\n",
                i, E, dos[i], meanr[i], sumr[i], countr[i]);
        }
                    
    fflush(f10);
    ENDDBG
    return 0;
    }
