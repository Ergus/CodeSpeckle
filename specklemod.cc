#include "specklemod.h"

speckle::speckle(int argc, char **argv):
    base_calculator(argc,argv),
    // Internal pointers very important
    thesolver(NULL),
    pointer_init(NULL),
    thehist(NULL),
    // ints
    start(1),                      
    end(1),    
    nov(4),
    Nscatterers(10),
    npmax(0),                      // This values are checked for errors
    save_interval(1),
    // double arrays
    VP(NULL),
    xpos(NULL),
 
    // doubles
    min(0.0), max(0.0),
    size(0.11),
    vDi(0.05635245901639344262),
    focal(40.0),
    vlambda(800.0e-6),
    vDi2(0.05635245901639344262),    
    focal2(30.0),    
    cofactor(1.0),    
    rphase(0.0),
    binsize(0.5),
    
    f18(NULL),    
    fprefix("output"),
    continuefile(""),

    // bools
    nrescale(false),    	   // shift and rescale sum2 speckle
    vectors(false),
    usegnuplot(false),    

    // complex number
    Vintensity(2.0*M_PI*M_PI),
    
    // internal arrays
    A(NULL), Vk(NULL),
    last_seed(0)

{
    STARTDBG;
    parse();

    // Some calculations after parsed the input.
    npxpu=npmax+1;
    Ntot=npmax*npmax*npmax;
    it=start;
    
    dx =size/double(npmax);
    dxi=double(npmax)/size;

    // This is the common name for the outputs, extension and
    // prefix can be added but this should not be modified.
    filename=fprefix+"_"+speckletype+"_"+timestr;
    
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

    if(!pointer_init){
        fprintf(stderr,"Error in assignement of init pointer\n");
        }
    
    // Some actions for rank 0 (or serial code have also rank==0)
    #ifdef DEBUG
    print();
    #endif

    // Construct the solver
    if(solvername=="mkl"){
        thesolver=new mkl_solver(Ntot,vectors,min,max,ncpu);
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
        thesolver=new plasma_solver(Ntot,vectors,min,max,ncpu,start_cpu);
        }
    #endif
    else{
        fprintf(stderr,"Solver=%s is not valid\n",solvername.c_str());
        fprintf(stderr,"Maybe a compilation time option not set\n");
        printme();
        exit(EXIT_FAILURE);
        }
    
    // Check for errors in the solver construction
    if (!thesolver){
        fprintf(stderr,"Error creating the solver class\n");
        printme();
        exit(EXIT_FAILURE);
        }

    ENDDBG;
    }

speckle::~speckle(){
    STARTDBG;
    if(thesolver) delete thesolver;
    if(thehist) delete thehist;
    if(xpos) free(xpos);    
    if(VP) free(VP);
    if(A) free(A);
    if(f18) fclose(f18);
    ENDDBG
    }

int speckle::init(int idseed){
    STARTDBG;
        
    //Allocate VP and xpos just ifthey are NULL  
    if(!VP) VP=(double*) malloc(npxpu*npxpu*npxpu*sizeof(double)); dbg_mem(VP);
    if(!xpos) xpos=(double*) malloc(npxpu*sizeof(double)); dbg_mem(xpos);
    
    //execute the init
    dbg((this->*pointer_init)(idseed));

    ENDDBG
    return(0);
    }

int speckle::ftspeckle(){
    STARTDBG
    // In all the cases N1==N2==N3 so the original Ni vars were sustituted by npx
    // The variables with a 2 attached means that are the original^2
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

    // The next for loop copies the data in the array VP to xreal, but without the extra columns for the boundary conditions
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
    double deltaT=(2.0*M_PI/size)*(2.0*M_PI/size);
    const int npx=npmax, npx2=npmax*npmax;    

    f18_printf("# The dimension of the matrix is %d\n", Ntot);
    
    int* vkx=(int*) malloc(Ntot*sizeof(int)); dbg_mem(vkx);
    int* vky=(int*) malloc(Ntot*sizeof(int)); dbg_mem(vky);
    int* vkz=(int*) malloc(Ntot*sizeof(int)); dbg_mem(vkz); 

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
    ENDDBG;
    return 0;
    }

int speckle::process_serial(int nvalues, double*array){
    STARTDBG;
    if(!thehist)
        thehist=new histogram(this);
    dbg(thehist->process(nvalues,array));
    ENDDBG;
    return 0;
    }


int speckle::calculate(int idseed){
    STARTDBG;
    last_seed=idseed;  //remember the last seed just to prevent errors
    
    #ifdef DEBUG  //This information is saved nly in debug mode
    char outputname[50];
    sprintf(outputname,"%s/%s_%d.out",
            dirname.c_str(),filename.c_str(),idseed);
    
    f18=fopen(outputname,"a"); dbg_mem(f18);
    print(f18,'#');
    f18_printf("# Seed\t %d\n",idseed);
    #endif // DEBUG
    
    dbg(init(idseed));
    //dbg(writespeckle());
    //dbg(correlationspeckle(idseed+1));
    dbg(ftspeckle());
    dbg(defineA());
    dbg(thesolver->solve(A));
    indices=thesolver->get_m();
    values=thesolver->get_w();
    #ifdef DEBUG
    fclose(f18); f18=NULL;
    #endif
    ENDDBG;
    return 0;
    }        

int speckle::correlationspeckle(int idseed){
    STARTDBG;
    srand(idseed); //Seed for random number generator
    
    const int Npart=150000,
              sizedouble=sizeof(double),
        sizeint=sizeof(int);
    
    int nthreads=1;  //this is a better approximation fot the best performance with omp    
    #ifdef _OPENMP
    nthreads=ncpu*0.8;    
    #endif // _OPENMP
    
    double *atau3    =(double*) calloc(nthreads*npmax,sizedouble); dbg_mem(atau3);
    double *atau3xy  =(double*) calloc(nthreads*npmax,sizedouble); dbg_mem(atau3xy);
    double *atau3z   =(double*) calloc(nthreads*npmax,sizedouble); dbg_mem(atau3z);
    double *atau3sq  =(double*) calloc(nthreads*npmax,sizedouble); dbg_mem(atau3sq);
    double *atau3xysq=(double*) calloc(nthreads*npmax,sizedouble); dbg_mem(atau3xysq);
    double *atau3zsq =(double*) calloc(nthreads*npmax,sizedouble); dbg_mem(atau3zsq);
    
    int *intau  =(int*) calloc(nthreads*npmax,sizeint); dbg_mem(intau);
    int *intauxy=(int*) calloc(nthreads*npmax,sizeint); dbg_mem(intauxy);
    int *intauz =(int*) calloc(nthreads*npmax,sizeint); dbg_mem(intauz);
    int *np     =(int*) calloc(3*Npart,sizeint); dbg_mem(np);
    
    const double size2=size*size, size2o4=0.25*size*size,
        stepcorfun = size/4000.0, vme=1.0;
    const int npx=npmax, npxo2=npmax/2, npxpu2=npxpu*npxpu;

    f18_printf("Computing Correlations\n");

    char namef13[512], namef19[512];
    sprintf(namef13,"%s/corr3d_%s_%d.dat",
            dirname.c_str(),filename.c_str(),idseed);
    sprintf(namef19,"%s/corfunanal3D_%s_%d.dat",
            dirname.c_str(),filename.c_str(),idseed);

    f18_printf("f13= %s\nf19= %s\n",namef13,namef19);

    FILE *f13=fopen(namef13,"w"); dbg_mem(f13);    
    FILE *f19=fopen(namef19,"w"); dbg_mem(f19);
    
    for(int i=0;i<Npart;i++){
        np[i*3  ]=(int)(frand()*size/dx);  //Here I didn't add 1 cause of the C indices
        np[i*3+1]=(int)(frand()*size/dx);
        np[i*3+2]=(int)(frand()*size/dx);
        }
    
    //i start in 1 because if i=j=0 nothing is made
    #pragma omp parallel num_threads(nthreads)
    {
    const int ithread = omp_get_thread_num();
    #ifdef DEBUG
    printf("Running correlation thread: %d process %d\n",ithread,rank);
    #endif // DEBUG
    #pragma omp for
    for(int i=0;i<Npart;i++){
        //double rx,ry,rz,rr,rxy;
        //int iabsiz;// ir;// iz;
        for(int j=0;j<i;j++){
            const double rx=dx*(np[i*3  ]-np[j*3  ]);
            const double ry=dx*(np[i*3+1]-np[j*3+1]);

            //is different because the int var is needed later
            const int iz=np[i*3+2]-np[j*3+2];   
            const double rz=dx*(iz);
            const double tmp1=VP[ np[j*3+2]*npxpu2+
                                  np[j*3+1]*npxpu +
                                  np[j*3  ]        ] - vme;

            const double rr=rx*rx+ry*ry+rz*rz;
            const double rxy=rx*rx+ry*ry;
            const int iabsiz=abs(iz);

            if(rr<size2o4){
                const int ir=(int)(sqrt(rr)*dxi);
                const double tmp2=VP[ np[i*3+2]*npxpu2+
                                      np[i*3+1]*npxpu +
                                      np[i*3  ]        ] - vme;                
                const double add = tmp1*tmp2;
                
                atau3[ithread*npmax+ir] += add;
                atau3sq[ithread*npmax+ir] += (add*add);
                intau[ithread*npmax+ir]++;
                }
            if(rxy<size2o4){
                const int ir=(int)(sqrt(rxy)*dxi);                
                const double tmp2=VP[ np[j*3+2]*npxpu2+
                                      np[i*3+1]*npxpu +
                                      np[i*3  ]        ] - vme;                
                const double add = tmp1*tmp2;                
                atau3xy[ithread*npmax+ir] += add;
                atau3xysq[ithread*npmax+ir] += (add*add);
                intauxy[ithread*npmax+ir]++;
                }
            if(iabsiz<npxo2){
                const int ir=(iz>0?iz:-iz);
                const double tmp2=VP[ np[i*3+2]*npxpu2+
                                      np[j*3+1]*npxpu +
                                      np[j*3  ]        ] - vme;
                
                const double add = tmp1*tmp2;                                
                atau3z[ithread*npmax+ir]+=add;
                atau3zsq[ithread*npmax+ir]+=(add*add);
                intauz[ithread*npmax+ir]++;
                }//end if
            }//end for
        }//end for
    //This is a reduction for the values before
    #ifdef _OPENMP
    #pragma omp for
    for(int i=0;i<npmax;i++){
        for(int j=1;j<nthreads;j++){
            atau3[i] += atau3[j*npmax+i] ;
            atau3sq[i] += atau3sq[j*npmax+i];
            intau[i] += intau[j*npmax+i] ;

            atau3xy[i] += atau3xy[j*npmax+i];
            atau3xysq[i] += atau3xysq[j*npmax+i];
            intauxy[i] += intauxy[j*npmax+i];

            atau3z[i]+=atau3z[j*npmax+i];
            atau3zsq[i]+=atau3zsq[j*npmax+i];
            intauz[i]+=intauz[j*npmax+i];     
            }
        }
    #endif // _OPENMP                     
        }// end openMP
    for(int i=0;i<=npxo2;i++){
        double val  =0.0, verr  =0.0, valxy=0.0, verrxy=0.0,
               valz =0.0, verrz =0.0, val2, vsq;
        if(intau[i] != 0){
            atau3[i]/=intau[i];
            val=atau3[i];
            val2=val*val;
            vsq=atau3sq[i]/intau[i];
            verr=(vsq-val2)/sqrt((double)(intau[i]));
            }
        if(intauxy[i] != 0){
            atau3xy[i]/=intauxy[i];
            valxy=atau3xy[i];
            val2=valxy*valxy;
            vsq=atau3xysq[i]/intauxy[i];
            verrxy=(vsq-val2)/sqrt((double)(intauxy[i]));
            }
        if(intauz[i] != 0){
            atau3z[i]/=intauz[i];
            valz=atau3z[i];
            val2=valz*valz;
            vsq=atau3zsq[i]/intauz[i];
            verrz=(vsq-val2)/sqrt((double)(intauz[i]));
            }
        fprintf(f13,"%lf %lf %lf %lf %lf %lf %lf\n",
                      ((double)i+0.5)*dx,val,verr,valxy,verrxy,valz,verrz);
        }//end for

    for(int i=1; i<=2000; i+=10){
        double xcorfun= stepcorfun*i,
            ltmp=(vDi*focal),
            argc   = M_PI*xcorfun/(vlambda*ltmp),
            c2circ  = 2.0*bessj1(argc)/(argc),
            fcorr3D= 3.0*(sin(argc)-argc*cos(argc))/(argc*argc*argc),
            argz = M_PI*xcorfun/(8.0*vlambda*ltmp*ltmp),
            czed=sin(argz)/argz;
        //all the previous vales are powered 2 in the original code
        //this is made in print time
        fprintf(f19,"%lf\t %lf\t %lf\t %lf\n",
                xcorfun, fcorr3D*fcorr3D, c2circ*c2circ, czed*czed );
        
        }//end for

    fclose(f13);
    fclose(f19);
    
    free(atau3);
    free(atau3xy);
    free(atau3z);
    free(atau3sq);
    free(atau3xysq);
    free(atau3zsq);
    
    free(intau);
    free(intauxy);
    free(intauz);
    free(np);
    ENDDBG;
    return 0;
    }


int speckle::writespeckle(){
    STARTDBG;

    if(!xpos || !VP){
        fprintf(stderr,"Error xpos and VP not initialized\n");
        printme();
        return -1;
        }
    
    char namef14[512];
    sprintf(namef14,"%s/speckle3d_%s_%d.dat",
            dirname.c_str(),filename.c_str(),last_seed);

    FILE *f14=fopen(namef14,"w"); dbg_mem(f14);
    const int npu=npxpu, npu2=npxpu*npxpu;
    
    for(int i=0;i<npxpu;i++){
        for(int j=0;j<npxpu;j++){
            for(int k=0;k<npxpu;k++){
                fprintf(f14,"%lf %lf %lf %lf\n",
                        xpos[k],xpos[j],xpos[i],VP[i*npu2+j*npu+k]);
                }
            fprintf(f14,"\n\n");
            }
        }
    fclose(f14);
    ENDDBG;
    return 0;
    }

inline double speckle::bessj1(double x){
    STARTDBG;
    double ax, z, xx,y,ans,ans1,ans2;
    // This is ugly but efficient.
    if((ax=fabs(x)) <8.0){
        y=x*x;
        ans1=x*(72362614232.0 +
                y*(-7895059235.0 +
                   y*(242396853.1 +
                      y*(-2972611.439 +
                         y*(15704.48260 +
                            y*(-30.16036606))))));
        ans2=144725228442.0+
            y*(2300535178.0+
               y*(18583304.74+
                  y*(99447.43394+
                     y*(376.9991397+
                        y*1.0))));
        ans=ans1/ans2;
        }
    else{
        z=8.0/ax;
        y=z*z;
        xx=ax-2.356194491;
        ans1=1.0+
            y*(0.183105e-2+
               y*(-0.3516396496e-4+
                  y*(0.2457520174e-5+
                     y*(-0.240337019e-6))));

        ans2=0.04687499995+
            y*(-0.2002690873e-3+
               y*(0.8449199096e-5+
                  y*(-0.88228987e-6+
                     y*0.105787412e-6)));
        if(x<0) ans=-sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
        else ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
        }
    ENDDBG;
    return ans;
    }
