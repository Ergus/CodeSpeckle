#include "specklemod.h"

speckle::speckle(int argc, char **argv):
    // Internal pointers very important
    thesolver(NULL),
    pointer_init(NULL),
    thehist(NULL),
    // ints
    nov(4),
    Nscatterers(10),
    npmax(0),                      // This values are checked for errors
    // double arrays
    VP(NULL),
    xpos(NULL),
 
    // doubles
    min(0), max(0),
    size(0.11),
    vDi(0.05635245901639344262),
    focal(40.0),
    vlambda(800.0e-6),
    vDi2(0.05635245901639344262),    
    focal2(30.0),    
    cofactor(1.0),    
    rphase(0.0),
    binsize(0.5),
    
    largc(argc), largv(argv),
    f18(NULL),    
    fprefix("NULL"),
    save_dir(""),
    continuefile(""),

    // bools
    nrescale(false),    	   // shift and rescale sum2 speckle
    vectors(false),
    usegnuplot(false),    

    // complex number
    Vintensity(2.0*M_PI*M_PI),
    
    // internal arrays
    A(NULL), Vk(NULL)
{
    STARTDBG
    parse();
    npxpu=npmax+1;
    Ntot=npmax*npmax*npmax;
    it=start;
    
    dx =size/double(npmax);
    dxi=double(npmax)/size;

    if (save_dir!="") {
        int created=mkdir(save_dir.c_str(),0777);
        if (created==0) printf("Process %d created dir %s\n",rank,save_dir.c_str());
        save_dir+="/";
        }
    if(fprefix=="NULL") fprefix="";   //prefix for output will beempty in this case
    
    //Set the speckle type to use
    if(speckletype=="spherical") pointer_init=&speckle::init_spherical;
    else if(speckletype=="sum2") pointer_init=&speckle::init_sum2;
    else if(speckletype=="shell") pointer_init=&speckle::init_shell;
    else if(speckletype=="single") pointer_init=&speckle::init_single;
    else{    
        fprintf(stderr,"Error the speckletype name provided is not valid\n");
        fprintf(stderr,"please check the input file or the command line\n");
        exit(EXIT_FAILURE);        
        }
    
    //set the solver to use
    if(solvername=="mkl"){
        thesolver=new mkl_solver(Ntot,vectors,min,max);
        }
    #ifdef UMAGMA
    else if(solvername=="magma"){
        thesolver=new magma_solver(Ntot,vectors,min,max,ngpu);
        }
    else if(solvername=="magma_2stage"){
        thesolver=new magma_solver_2stage(Ntot,vectors,min,max,ngpu);
        }
    #endif
    #ifdef UPLASMA
    else if(solvername=="plasma"){
        thesolver=new plasma_solver(Ntot,vectors,min,max);
        }
    #endif
    else{
        fprintf(stderr,"Solver=%s is not valid\n",solvername.c_str());
        fprintf(stderr,"Maybe a compilation time option not set\n");
        printme();
        exit(EXIT_FAILURE);
        }
    #ifdef DEBUG
    print();
    #endif //DEBUG
    ENDDBG
    }

speckle::~speckle(){
    STARTDBG
    if(thesolver) delete thesolver;
    if(thehist) delete thehist;
    if(xpos) free(xpos);    
    if(VP) free(VP);
    if(A) free(A);
    if(f18) fclose(f18);
    ENDDBG
    }


int speckle::init(int idseed){
    STARTDBG
    //Allocate VP and xpos just ifthey are NULL  
    if(!VP) dbg_mem((VP=(double*) malloc(npxpu*npxpu*npxpu*sizeof(double))));
    if(!xpos) dbg_mem((xpos=(double*) malloc(npxpu*sizeof(double))));

    #ifdef DEBUG
    if(fprefix!=""){
        sprintf(outputname,"%s%s%s_%d.out",
                save_dir.c_str(),fprefix.c_str(),speckletype.c_str(),idseed);
    
        f18=fopen(outputname,"a"); dbg_mem(f18);
        print(f18,'#');
        fprintf(f18,"# Seed\t %d\n",idseed);
        }        
            
    #endif // DEBUG

    //execute the init
    (this->*pointer_init)(idseed);
    ENDDBG
    return(0);
    }

int speckle::ftspeckle(){
    STARTDBG
    //In all the cases N1==N2==N3 so the original Ni vars were sustituted by npx
    const int npx=npmax, npx2=npx*npx, npu=npxpu, npu2=npxpu*npxpu, padx=npx/2+1;
    
    MKL_LONG N[]={ npx, npx, npx};
    //MKL_LONG rstrides[]={ 0,     npx2,  npx, 1};
    //MKL_LONG cstrides[]={ 0, npx*padx, padx, 1};
    MKL_LONG rstrides[]={ 0,      npx, npx2, 1};
    MKL_LONG cstrides[]={ 0, npx*padx, padx, 1};

    //This 2 arrays will have dimension 3
    f18_printf("# Allocate data arrays x_real\n");
    double *x_real=(double*) malloc(npx*npx2*sizeof(double)); dbg_mem(x_real);
    double complex *x_cmplx=(double complex*) malloc(padx*npx2*sizeof(double complex));
    dbg_mem(x_cmplx);
    
    if(Vk) free(Vk);
    Vk=(double complex *) malloc(npx*npx2*sizeof(double complex)); dbg_mem(Vk);

    //The next for loop copies the data in the array VP to xreal, but without the extra columns for the boundary conditions
    f18_printf("# Initialize data for real-to-complex FFT\n");
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
    
    f18_printf("# Create DFTI descriptor for real transform\n");
    dbg(DftiCreateDescriptor(&hand,DFTI_DOUBLE,DFTI_REAL,3,N));

    f18_printf("# Set out-of-place\n");
    dbg(DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE));

    f18_printf("# Set CCE storage\n");
    dbg(DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX));
    dbg(DftiSetValue(hand, DFTI_FORWARD_SCALE, scaleforward));
    dbg(DftiSetValue(hand, DFTI_INPUT_STRIDES, rstrides));
    dbg(DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cstrides));

    f18_printf("# Commit DFTI descriptor\n");
    dbg(DftiCommitDescriptor(hand));

    f18_printf("# Compute forward transform\n");
    dbg(DftiComputeForward(hand, x_real, x_cmplx));
                   
    f18_printf("# Reconfigure DFTI descript for back transform\n");
    dbg(DftiSetValue(hand, DFTI_INPUT_STRIDES, cstrides));
    dbg(DftiSetValue(hand, DFTI_OUTPUT_STRIDES, rstrides));
    dbg(DftiSetValue(hand,DFTI_BACKWARD_SCALE, scalebackward));
    dbg(DftiCommitDescriptor(hand));
    dbg(DftiComputeBackward(hand, x_cmplx, x_real));

    f18_printf("# Verify the result after a forward and backward Fourier Transforms\n");
    int count=0;
    for(int i=0;i<npx;i++){
        for(int j=0;j<npx;j++){
            for(int k=0;k<npx;k++){
                count+=(fabs(x_real[i*npx2+j*npx+k]-VP[i*npu2+j*npu+k])>0.01 ? 1 : 0);
                }
            }
        }

    if(count>0){
        fprintf(stderr,"WARNING!!!======\n");
        fprintf(stderr,"Elements with an error > 0.01: %d\n",count);
        }


    //    show(x_cmplx,);
    //This is completly modified way to do this because no mod or extra if are needed
    for(int i=0; i<npx; i++){
        int ni=(npx-i)%npx;
        for(int j=0; j<npx; j++){
            int nj=(npx-j)%npx;
            for(int k=0;k<padx;k++){ //same vales than padded
                Vk[i*npx2+j*npx+k]=Vintensity*x_cmplx[i*npx*padx+j*padx+k];
                }
            for(int k=padx; k<npx; k++){
                int nk=(npx-k)%npx;
                Vk[i*npx2+j*npx+k]=Vintensity*conj(x_cmplx[ni*npx*padx+nj*padx+nk]);
                }
            }
        }
    
    free(x_real);
    free(x_cmplx);
    f18_printf("# Ftspeckle ended ok\n");
    ENDDBG
    return 0;
    }

int speckle::defineA(){
    STARTDBG
    double deltaT=(2.0*M_PI/size);
    const int npx=npmax, npx2=npmax*npmax;    

    f18_printf("# The dimension of the matrix is %d\n", Ntot);
    
    int* vkx=(int*) malloc(Ntot*sizeof(int)),
       * vky=(int*) malloc(Ntot*sizeof(int)),
       * vkz=(int*) malloc(Ntot*sizeof(int));
    //If error in this line check the 3 lines before
    dbg_mem(vkx); dbg_mem(vky); dbg_mem(vkz); 

    double* diagTk=(double*) malloc(Ntot*sizeof(double)); dbg_mem(diagTk);
    
    if(A) free(A);
    A=(double complex *) calloc(Ntot*Ntot,sizeof(double complex)); dbg_mem(A);
    
    f18_printf("# The dimension of the matrix is %d\n", Ntot);

    //----------Here we define the matrix A-----------------
    //vectors of the k-components and of the kinetic energy (in units of deltaT=hbar^2*unitk^2/2m)
    int ntk=0;
    int kz=-npmax/2;
    for(int nkz=0;nkz<npmax;nkz++){
        int ky=-npmax/2;
        for(int nky=0;nky<npmax; nky++){
            int kx=-npmax/2;
            for(int nkx=0 ;nkx<npmax; nkx++){
                diagTk[ntk]=kx*kx+ky*ky+kz*kz;
                vkx[ntk]=kx;
                vky[ntk]=ky;
                vkz[ntk]=kz;
                ntk++;
                kx++;
                }
            ky++;
            }
        kz++;
        }
        
    //This loop is transversed access, but this is the way
    //it is made in the original code.
    //This can be made much more efficiently in the previous loop, but
    //maybe this way is usefull or readable.
    f18_printf("# Definition of A\n");
    for(int i=0;i<Ntot;i++){
        for(int j=i+1;j<Ntot;j++){
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
    free(Vk); Vk=NULL;//Vk is not used anymore and we need memory to diagonalization
    f18_printf("# Matrix A Defined correctly\n");
    ENDDBG
    return 0;
    }


int speckle::process(int nvalues, double*array){
    STARTDBG
    if(!thehist)
        thehist=new histogram(this); dbg_mem(thehist);
    dbg(thehist->process(nvalues,array));
    ENDDBG
    return 0;
    }

