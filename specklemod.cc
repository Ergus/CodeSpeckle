#include "specklemod.h"
#include "speckleauxiliary.h"

speckle::speckle(int argc, char** argv):
    //Secound constructor using standard input.
    Nscatterers(100),
    nrescale(false),    	  //shift and rescale sum2 speckle
    size(0.11),
    vDi(0.05635245901639344262),
    focal(40.0),
    vlambda(800.0e-6),
    vDi2(0.05635245901639344262),
    focal2(30.0),
    cofactor(1.0),
    rphase(0.0),
    npmax(0),                 //This values are checked for errors
    nov(4),
    Vk(NULL),
    A(NULL),
    largc(argc),
    largv(argv){
    
    struct option localopts[] = {
        {"Nscatterers", required_argument, 0, 'N'},
        {"size", required_argument, 0, 's'},
        {"vDi", required_argument, 0, 'v'},
        {"vDi2", required_argument, 0, 'V'},
        {"focal", required_argument, 0, 'f'},
        {"focal2", required_argument, 0, 'F'},
        {"vlambda", required_argument, 0, 'l'},
        {"cofactor", required_argument, 0, 'c'},
        {"rphase", required_argument, 0, 'p'},
        {"npmax", required_argument, 0, 'm'},
        {"speckletype", required_argument, 0, 't'},
        {"nov", required_argument, 0, 'n'},
                
        {"rescale", no_argument, 0, 'r'},
        {"help", no_argument, 0, 'h'},
                
        {0, 0, 0, 0}
        };

    longopts=localopts;
    parser();
    npxpu=npmax+1;    
    
    VP=(double*) malloc(npxpu*npxpu*npxpu*sizeof(double));
    dbg_mem(VP);
    xpos=(double*) malloc(npxpu*sizeof(double));
    dbg_mem(xpos);

    dx =size/double(npmax);
    dxi=double(npmax)/size;

    print();
    }

int speckle::init(int idseed){
    switch(speckletype){
        case spherical:
            init_spherical(idseed);
            break;
        case sum2:
            init_sum2(idseed);
            break;
        case single:
            init_single(idseed);
            break;
        case shell:
            init_shell(idseed);
            break;
        default:
            fprintf(stderr,"Error: calling %s in %s:%d\n",
                    __PRETTY_FUNCTION__, __FILE__,__LINE__);
            
        }
    }

speckle::~speckle(){
    if(xpos) free(xpos);    
    if(VP) free(VP);
    if(A) free(A);
    }

double speckle::EXVT(double xco,double yco,double zco){
    double ya[nov*nov*nov], y, dy;
    int ix=(int)(xco * dxi),
        iy=(int)(yco * dxi),
        iz=(int)(zco * dxi),
        last=npxpu-nov,
        half=(nov-1)/2,
        nov2=nov*nov,
        npxpu2=npxpu*npxpu,
        kx=min(max(ix-half,0),last),
        ky=min(max(iy-half,0),last),
        kz=min(max(iz-half,0),last);
    
    for(int i=kz; i<kz+nov; i++){
        for(int j=ky; j<ky+nov; j++){
            for(int k=kx; k<kx+nov; k++){
                ya[ (i-kz)*nov2
                   +(j-ky)*nov
                   +(k-kx)] = VP[ i*npxpu2
                                 +j*npxpu
                                 +k];
                }
            }
        }

    polin3(&xpos[kx],&xpos[ky],&xpos[kz],
           ya, nov, nov, nov,xco,yco,zco,y,dy);
        
    return y;
    }

void speckle::correlationspeckle(int idseed){
    const int Npart=150000,
              npartmax=200000,
              sizedouble=sizeof(double),
              sizeint=sizeof(int);
    
    double *atau3    =(double*) calloc(npmax,sizedouble),
           *atau3xy  =(double*) calloc(npmax,sizedouble),
           *atau3z   =(double*) calloc(npmax,sizedouble),
           *atau3sq  =(double*) calloc(npmax,sizedouble),
           *atau3xysq=(double*) calloc(npmax,sizedouble),
           *atau3zsq =(double*) calloc(npmax,sizedouble);
    
    int *intau  =(int*) calloc(npmax,sizeint),
        *intauxy=(int*) calloc(npmax,sizeint),
        *intauz =(int*) calloc(npmax,sizeint),
        *np     =(int*) calloc(3*npartmax,sizeint);
    
    double vme=1.0;
    
    const double size2=size*size, size2o4=0.25*size*size, stepcorfun = size/4000.;
    const int npx=npmax, npxo2=npmax/2, npxpu2=npxpu*npxpu;
    
    int nx,usedsize;

    printf("Computing Correlations\n");

    FILE *f13=fopen("specklecorr3D.dat","a");
    FILE *f19=fopen("specklecorfunanal3D.dat","a");

    for(int i=0;i<Npart;i++){
        np[i*3  ]=(int)(frand()*size/dx);  //Here I didn't add 1 cause of the C indices
        np[i*3+1]=(int)(frand()*size/dx);
        np[i*3+2]=(int)(frand()*size/dx);
        }
    
    //i start in 1 because if i=j=0 nothing is made
    for(int i=1;i<Npart;i++){
        double rx,ry,rz,rr,rxy;
        int iabsiz, ir, iz;
        for(int j=0;j<i;j++){
            rx=dx*(np[i*3  ]-np[j*3  ]);
            ry=dx*(np[i*3+1]-np[j*3+1]);
            
            iz=np[i*3+2]-np[j*3+2];   //is different because the int var is needed later
            rz=dx*(iz);
            
            const double tmp=VP[np[j*3+2]*npxpu2 + np[j*3+1]*npxpu + np[j*3]] - vme,
                add = tmp*tmp,
                add2= add*add;
            
            rr=rx*rx+ry*ry+rz*rz;
            rxy=rx*rx+ry*ry;
            iabsiz=abs(iz);

            if(rr<size2o4){
                ir=(int)(sqrt(rr)/dx);
                atau3[ir] += add;
                atau3sq[ir] += add2;
                intau[ir] += + 1;
                }
            if(rxy<size2o4){
                ir=(int)(sqrt(rxy)/dx);
                atau3xy[ir] += add;
                atau3xysq[ir] += add2;
                intauxy[ir]++;
                }
            if(iabsiz<npxo2){
                ir=iabsiz;
                atau3z[ir]+=add;
                atau3zsq[ir]+=add2;
                intauz[ir]++;
                }//end if
            }//end for
        }//end for

    for(int i=0;i<=npxo2;i++){
        double val  =0.0, verr  =0.0, valxy=0.0, verrxy=0.0,
               valz =0.0, verrz =0.0, val2, vsq;
        if(intau[i] != 0){
            val=(atau3[i]/=intau[i]);
            val2=val*val;
            vsq=atau3sq[i]/intau[i];
            verr=(vsq-val2)/sqrt((double)(intau[i]));
            }
        if(intauxy[i] != 0){
            valxy=(atau3xy[i]/=intauxy[i]);
            val2=val*val;
            vsq=atau3xysq[i]/intauxy[i];
            verrxy=(vsq-val2)/sqrt((double)(intauxy[i]));
            }
        if(intauz[i] != 0){
            valz=(atau3z[i]/=intauz[i]);
            val2=val*val;
            vsq=atau3zsq[i]/intauz[i];
            verrz=(vsq-val2)/sqrt((double)(intauz[i]));
            }
        fprintf(f13,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",
                      ((double)i+0.5)*dx,val,verr,valxy,verrxy,valz,verrz);
        }//end for

    for(int i=1; i<=2000; i+=10){
        double xcorfun= stepcorfun*i,
               argc   = M_PI*xcorfun/(vDi*vlambda*focal),
               c2circ  = 2.0*bessj1(argc)/(argc),
               fcorr3D= 3.0*(sin(argc)-argc*cos(argc))/(argc*argc*argc),
               ltmp=(vDi*focal),
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
    }

int speckle::ftspeckle(){
    //In all the cases N1==N2==N3 so the original Ni vars in the code were sustituted for npx
    const int npx=npmax, npx2=npx*npx, npu=npxpu, npu2=npxpu*npxpu, padx=npx/2+1;
    
    MKL_LONG status,
             N[3]={ npx, npx, npx},
             rstrides[4]={ 0,     npx2,  npx, 1},
             cstrides[4]={ 0, npx*padx, padx, 1};

    //This 2 arrays will have dimension 3
    printf("Allocate data arrays\n");             
    double *x_real=(double*) malloc(npx*npx2*sizeof(double)); dbg_mem(x_real);
    double complex *x_cmplx=(double complex*) malloc(padx*npx2*sizeof(double complex));
    dbg_mem(x_cmplx);
    
    if(Vk) free(Vk);
    Vk=(double complex *) malloc(npx*npx2*sizeof(double complex)); dbg_mem(Vk);

    //The next for loop copies the data in the array VP to xreal, but without the extra columns for the boundary conditions
    printf("Initialize data for real-to-complex FFT\n");
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
    
    printf("Create DFTI descriptor for real transform\n");
    dbg(DftiCreateDescriptor(&hand,DFTI_DOUBLE,DFTI_REAL,3,N));

    printf("Set out-of-place\n");
    dbg(DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE));

    printf("Set CCE storage\n");
    dbg(DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX));

    dbg(DftiSetValue(hand, DFTI_FORWARD_SCALE, scaleforward));

    dbg(DftiSetValue(hand, DFTI_INPUT_STRIDES, rstrides));

    dbg(DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cstrides));
 
    printf("Commit DFTI descriptor\n");
    dbg(DftiCommitDescriptor(hand));

    printf("Compute forward transform\n");
    dbg(DftiComputeForward(hand, x_real, x_cmplx));
                   
    printf("Reconfigure DFTI descript for back transform\n");
    dbg(DftiSetValue(hand, DFTI_INPUT_STRIDES, cstrides));

    dbg(DftiSetValue(hand, DFTI_OUTPUT_STRIDES, rstrides));

    dbg(DftiSetValue(hand,DFTI_BACKWARD_SCALE, scalebackward));

    dbg(DftiCommitDescriptor(hand));

    dbg(DftiComputeBackward(hand, x_cmplx, x_real));

    printf("Verify the result after a forward and backward Fourier Transforms\n");
    int count=0;
    for(int i=0;i<npx;i++){
        for(int j=0;j<npx;j++){
            for(int k=0;k<npx;k++){
                count+=(fabs(x_real[i*npx2+j*npx+k]-VP[i*npu2+j*npu+k])>0.01);
                }
            }
        }

    if(count>0){
        printf(" WARNING!!!======\n");
        printf(" Elements with an error > 0.01: %d\n",count);
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
    printf("ftspeckle ended ok\n");
    return 0;
    }

void speckle::defineA(){
    double deltaT=(2.0*M_PI/size);
    const int Ntot=npmax*npmax*npmax, npx=npmax, npx2=npmax*npmax;    

    printf("The dimension of the matrix is %d\n", Ntot);
    
    int* vkx=(int*) malloc(Ntot*sizeof(double)),
       * vky=(int*) malloc(Ntot*sizeof(double)),
       * vkz=(int*) malloc(Ntot*sizeof(double));

    double* diagTk=(double*) malloc(Ntot*sizeof(double));
    
    if(A) free(A);
    A=(double complex *) calloc(Ntot*Ntot,sizeof(double complex)); dbg_mem(A);
    
    printf("The dimension of the matrix is %d\n", Ntot);

    //----------Here we define the matrix A-----------------
    //vectors of the k-components and of the kinetic energy (in units of deltaT=hbar^2*unitk^2/2m)    
    for(int nkx=0,kx=-npmax/2,ntk=0 ;nkx<npmax; nkx++,kx++){         //ntk is declared here
        for(int nky=0,ky=-npmax/2 ;nky<npmax; nky++,ky++){
            for(int nkz=0,kz=-npmax/2; nkz<npmax; nkz++,kz++,ntk++){ //but incremented here
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
    printf("Definition of A\n");    
    for(int i=0;i<Ntot;i++){
        for(int j=i+i;j<Ntot;j++){
            int kdiffx=(vkx[j]-vkx[i]+npmax)%npmax;
            int kdiffy=(vky[j]-vky[i]+npmax)%npmax;
            int kdiffz=(vkz[j]-vkz[i]+npmax)%npmax;
            double complex t=Vk[kdiffz*npx2+kdiffy*npx+kdiffx];
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
    printf("Matrix A Defined correctly\n");
    }
