#include "specklemod.h"


speckle::speckle(parser *thepar):
    Vk(NULL),
    xpos(NULL),
    A(NULL),
    VP(NULL),
    //now copy the values needed from the parser
    Nscatterers(thepar->Nscatterers),
    npmax(thepar->npmax),
    nov(thepar->nov),
    size(thepar->size),
    vDi(thepar->vDi),
    vDi2(thepar->vDi2),    
    focal(thepar->focal),
    focal2(thepar->focal2),    
    vlambda(thepar->vlambda),
    cofactor(thepar->cofactor),
    rphase(thepar->rphase),
    nrescale(thepar->nrescale),
    vectors(thepar->vectors),
    pointer_init(NULL),
    parptr(thepar),   //The parser pointer will point to the argument parser
    f18(NULL),
    Vintensity(thepar->Vintensity),
    fprefix(thepar->fprefix),
    speckletype(thepar->speckletype),
    solver(thepar->solver),
    min(thepar->min),
    max(thepar->max),
    thesolver(NULL)
{
    
    npxpu=npmax+1;
    Ntot=npmax*npmax*npmax;

    dx =size/double(npmax);
    dxi=double(npmax)/size;

    sprintf(outputname,"/dev/null");    

    //Set the speckle type to use
    if(speckletype=="spherical") pointer_init=&speckle::init_spherical;
    else if(speckletype=="sum2") pointer_init=&speckle::init_sum2;
    else if(speckletype=="shell") pointer_init=&speckle::init_shell;
    else if(speckletype=="single") pointer_init=&speckle::init_single;
    else{
        if(parptr->rank==0){
            fprintf(stderr,"Error the speckletype name provided is not valid\n");
            fprintf(stderr,"please check the input file or the command line\n");
            }
        exit(EXIT_FAILURE);
        }
    
    //set the solver to use
    if(solver=="mkl"){
        thesolver=new mkl_solver(Ntot,vectors,min,max);
        }
    #ifdef UMAGMA
    else if(solver=="magma"){
        thesolver=new magma_solver(Ntot,vectors,min,max,2);
        }
    else if(solver=="magma_2stage"){
        thesolver=new magma_solver_2stage(Ntot,vectors,min,max,2);
        }
    #endif
    #ifdef UPLASMA
    else if(solver=="plasma"){
        thesolver=new plasma_solver(Ntot,vectors,min,max);
        }
    #endif
    else{
        fprintf(stderr,"The value provided for the solver=%s is not valid\n",solver.c_str());
        printme();
        exit(EXIT_FAILURE);
        }
    }

int speckle::init(int idseed){
    //Allocate VP and xpos just ifthey are NULL  
    if(!VP) dbg_mem((VP=(double*) malloc(npxpu*npxpu*npxpu*sizeof(double))));
    if(!xpos) dbg_mem((xpos=(double*) malloc(npxpu*sizeof(double))));

    //set the output name for this seed
    if((fprefix!="null")&&(fprefix!="NULL")){
        sprintf(outputname,"%s_%s_%d.out",fprefix.c_str(),speckletype.c_str(),idseed);
        }
    //open the file just for a moment the close it
    f18=fopen(outputname,"a"); dbg_mem(f18);
    parptr->print(f18,'#');
    fprintf(f18,"# seed\t %d\n",idseed);
    //execute the init
    (this->*pointer_init)(idseed);
    fclose(f18);
    return(0);
    }

speckle::~speckle(){
    if(thesolver) delete thesolver;
    if(xpos) free(xpos);    
    if(VP) free(VP);
    if(A) free(A);
    }


int speckle::ftspeckle(){
    //In all the cases N1==N2==N3 so the original Ni vars in the code were sustituted for npx
    const int npx=npmax, npx2=npx*npx, npu=npxpu, npu2=npxpu*npxpu, padx=npx/2+1;
    
    MKL_LONG status;
    MKL_LONG N[]={ npx, npx, npx};
    MKL_LONG rstrides[]={ 0,     npx2,  npx, 1};
    MKL_LONG cstrides[]={ 0, npx*padx, padx, 1};

    //This 2 arrays will have dimension 3
    f18=fopen(outputname,"a");
    fprintf(f18,"# Allocate data arrays x_real\n");
    double *x_real=(double*) malloc(npx*npx2*sizeof(double)); dbg_mem(x_real);
    double complex *x_cmplx=(double complex*) malloc(padx*npx2*sizeof(double complex));
    dbg_mem(x_cmplx);
    
    if(Vk) free(Vk);
    Vk=(double complex *) malloc(npx*npx2*sizeof(double complex)); dbg_mem(Vk);

    //The next for loop copies the data in the array VP to xreal, but without the extra columns for the boundary conditions
    fprintf(f18,"# Initialize data for real-to-complex FFT\n");
    for(int i=0;i<npx;i++){
        for(int j=0;j<npx;j++){
            for(int k=0;k<npx;k++){
                x_real[i*npx2+j*npx+k]=VP[i*npu2+j*npu+k];
                }
            }
        }

	DFTI_DESCRIPTOR_HANDLE hand=0;
    double scaleforward=1.0/(npx*npx2),
           scalebackward=1.0;
    
    fprintf(f18,"# Create DFTI descriptor for real transform\n");
    dbg(DftiCreateDescriptor(&hand,DFTI_DOUBLE,DFTI_REAL,3,N));

    fprintf(f18,"# Set out-of-place\n");
    dbg(DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE));

    fprintf(f18,"# Set CCE storage\n");
    dbg(DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX));

    dbg(DftiSetValue(hand, DFTI_FORWARD_SCALE, scaleforward));

    dbg(DftiSetValue(hand, DFTI_INPUT_STRIDES, rstrides));

    dbg(DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cstrides));
 
    fprintf(f18,"# Commit DFTI descriptor\n");
    dbg(DftiCommitDescriptor(hand));

    fprintf(f18,"# Compute forward transform\n");
    dbg(DftiComputeForward(hand, x_real, x_cmplx));
                   
    fprintf(f18,"# Reconfigure DFTI descript for back transform\n");
    dbg(DftiSetValue(hand, DFTI_INPUT_STRIDES, cstrides));

    dbg(DftiSetValue(hand, DFTI_OUTPUT_STRIDES, rstrides));

    dbg(DftiSetValue(hand,DFTI_BACKWARD_SCALE, scalebackward));

    dbg(DftiCommitDescriptor(hand));

    dbg(DftiComputeBackward(hand, x_cmplx, x_real));

    fprintf(f18,"# Verify the result after a forward and backward Fourier Transforms\n");
    int count=0;
    for(int i=0;i<npx;i++){
        for(int j=0;j<npx;j++){
            for(int k=0;k<npx;k++){
                count+=(fabs(x_real[i*npx2+j*npx+k]-VP[i*npu2+j*npu+k])>0.01);
                }
            }
        }

    if(count>0){
        fprintf(stderr,"WARNING!!!======\n");
        fprintf(stderr,"Elements with an error > 0.01: %d\n",count);
        }
    
    //This is completly modified way to do this because no mod or extra if are needed
    for(int i=0, ni=npx; i<npx; i++, ni--){
        for(int j=0, nj=npx; j<npx; j++, nj--){
            for(int k=0;k<padx;k++){ //same vales than padded
                Vk[i*npx2+j*npx+k]=Vintensity*x_cmplx[i*npx2+j*npx+k];
                }
            for(int k=padx, nk=npx-padx; k<npx; k++, nk--){
                Vk[i*npx2+j*npx+k]=Vintensity*conj(x_cmplx[ni*npx*padx+nj*padx+nk]);
                }
            }
        }
    
    free(x_real);
    free(x_cmplx);
    fprintf(f18,"# Ftspeckle ended ok\n");
    fclose(f18);
    return 0;
    }

int speckle::defineA(){
    double deltaT=(2.0*M_PI/size);
    const int npx=npmax, npx2=npmax*npmax;    

    f18=fopen(outputname, "a");
    fprintf(f18,"# The dimension of the matrix is %d\n", Ntot);
    
    int* vkx=(int*) malloc(Ntot*sizeof(double)),
       * vky=(int*) malloc(Ntot*sizeof(double)),
       * vkz=(int*) malloc(Ntot*sizeof(double));

    double* diagTk=(double*) malloc(Ntot*sizeof(double)); dbg_mem(diagTk);
    
    if(A) free(A);
    A=(double complex *) calloc(Ntot*Ntot,sizeof(double complex)); dbg_mem(A);
    
    fprintf(f18,"# The dimension of the matrix is %d\n", Ntot);

    //----------Here we define the matrix A-----------------
    //vectors of the k-components and of the kinetic energy (in units of deltaT=hbar^2*unitk^2/2m)    
    for(int nkx=0,kx=-npmax/2,ntk=0 ;nkx<npmax; nkx++,kx++){  //ntk is declared here
        for(int nky=0,ky=-npmax/2 ;nky<npmax; nky++,ky++){
            for(int nkz=0,kz=-npmax/2; nkz<npmax; nkz++,kz++,ntk++){ //but increment here
                diagTk[ntk]=kx*kx+ky*ky+kz*kz;
                vkx[ntk]=kx;
                vky[ntk]=ky;
                vkz[ntk]=kz;
                }
            }
        }

    //This loop is transversed access, but this is the way
    //it is made in the original code.
    //This can be made much more efficiently in the previous loop, but
    //maybe this way is usefull or readable.
    fprintf(f18,"# Definition of A\n");    
    for(int i=0;i<Ntot;i++){
        for(int j=i+i;j<Ntot;j++){
            int kdiffx=(vkx[j]-vkx[i]+npmax)%npmax;
            int kdiffy=(vky[j]-vky[i]+npmax)%npmax;
            int kdiffz=(vkz[j]-vkz[i]+npmax)%npmax;
            const double complex t=Vk[kdiffz*npx2+kdiffy*npx+kdiffx];
            A[j*Ntot+i]=t;
            A[i*Ntot+j]=conj(t);
            }
        A[i*Ntot+i]=Vk[0]+(deltaT*diagTk[i]);
        }
    
    free(diagTk);
    free(vkx);
    free(vky);
    free(vkz);
    free(Vk); //Vk is not used anymore and we need memory to diagonalization
    fprintf(f18,"# Matrix A Defined correctly\n");
    fclose(f18);
    return 0;
    }
