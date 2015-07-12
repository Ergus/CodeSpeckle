#include "specklemod.h"
#include "speckleauxiliary.h"

char Specklename[][10]={"spherical","sum2","single","shell"};

speckle::speckle(int onpmax, Speckletype ospeckletype):
    // variables to define the speckle 
    Nscatterers(100),
    nrescale(false),    	//shift and rescale sum2 speckle
    size(0.11),
    vDi(0.05635245901639344262),
    focal(40.0),
    vlambda(800.0e-6),
    vDi2(0.05635245901639344262),
    focal2(30.0),
    cofactor(1.0),
    rphase(0.0),
    VP(NULL),
    xpos(NULL),
    npmax(onpmax),
    npxpu(onpmax+1),
    speckletype(ospeckletype),
    nov(4),
    Vk(NULL){
        VP=(double*) malloc(npxpu*npxpu*npxpu*sizeof(double));
        xpos=(double*) malloc(npxpu*sizeof(double));
        x=(double*) malloc(npmax*sizeof(double));

        dx =size/double(npmax);
        dxi=double(npmax)/size;
        }

int speckle::init(int idseed){
    //For to be used internally, not needed but for to make
    //the code looks like the Fortran One
    const int npx=npmax, npx2=npx*npx, npu=npxpu, npu2=npxpu*npxpu;
    
    const int ngridvI =1000;
    double *pvI=(double*) calloc(ngridvI,sizeof(double));

    double vmaxradscat =1.0/(2*vDi ),
           vmaxradscat2=1.0/(2*vDi2),
           deltaR =vlambda*focal /size,
           deltaR2=vlambda*focal2/size,        
           //for to avoid this extra multiplication too many times
           vmaxradscatquad=vmaxradscat*vmaxradscat,
           vmaxradscat2quad=vmaxradscat2*vmaxradscat2,
           aver=0.0,
           invlambdafoc=1.0/(vlambda*focal),
           invlambdafoc2,               //this is only references in plain and sum2
           stdev=0.0,
           vImax=10.0,
           dvI=vImax/ngridvI,
           adpp=1.0/(npx2*npx),
           xcord[3];
    
    double *xscat=(double*) malloc(Nscatterers*sizeof(double)),
           *yscat=(double*) malloc(Nscatterers*sizeof(double)),
           *alphascat=(double*) malloc(Nscatterers*sizeof(double)),
           *betascat=(double*) malloc(Nscatterers*sizeof(double)),
           *zscat=NULL,               //spherical and plane
           *xscat2=NULL,              //for sum
           *zscat2=NULL,              //for sum
           *alphascat2=NULL,          //for sum
           *betascat2=NULL;           //for sum
    
    double complex ffaux;             //this is used as a temp, but in different loops
    //This will have dimension 3
    double complex *ff3=(double complex*) calloc(npmax*npmax*npmax,sizeof(double complex)),
                   *ff3_2=NULL;

    printf("Creating Speckle\n");
    FILE *f18=fopen("specklepi3D.dat","a");
    srand(idseed); //Seed for random number generator

    for(int i=0;i<npx;i++){
        double tmp=i*dx;
        x[i]=tmp;
        xpos[i]=tmp;
        }
    xpos[npx]=size;
    printf("N grid = %d\n",npx);

    //sample moduli and phases of scattered waves
    for(int i=0;i<Nscatterers;i++){
        alphascat[i] = 2.0*frand();
        betascat[i] = 2.0*M_PI*frand()-M_PI;
        }

    if(speckletype==sum2){
        xscat2=(double*) malloc(Nscatterers*sizeof(double));
        zscat2=(double*) malloc(Nscatterers*sizeof(double));
        alphascat2=(double*) malloc(Nscatterers*sizeof(double));
        betascat2=(double*) malloc(Nscatterers*sizeof(double));
        ff3_2=(double complex*) calloc(npmax*npmax*npmax,sizeof(double complex));
        for(int i=0;i<Nscatterers;i++){
            alphascat2[i] = 2.0*frand();
            betascat2[i] = 2.0*M_PI*frand()-M_PI;
            }
        }
    
    //sample positions of random scatterers
    if( (speckletype==single) || (speckletype==sum2) ){
        printf("Speckletype: %s\n",Specklename[speckletype]);
        double tmp1, tmp2;
        for(int i=0;i<Nscatterers;i++){
            do{
                tmp1 = 2.0*vmaxradscat*frand()-vmaxradscat;
                tmp2 = 2.0*vmaxradscat*frand()-vmaxradscat;
                //putting points in grid
                //Here round remplaces ANINT in fortran code
                xscat[i] = deltaR * round(tmp1/deltaR);
                yscat[i] = deltaR * round(tmp2/deltaR);
                }while((xscat[i]*xscat[i]+yscat[i]*yscat[i])>vmaxradscatquad);
            }
        //This is only for the sum2 case, almost the same code above
        if(speckletype==sum2){
            for(int i=0;i<Nscatterers;i++){
                do{
                    tmp1 = 2.0*vmaxradscat2*frand()-vmaxradscat2;
                    tmp2 = 2.0*vmaxradscat2*frand()-vmaxradscat2;
                    //putting points in grid
                    //Here round remplaces ANINT in fortran code
                    xscat2[i] = deltaR2 * round(tmp1/deltaR);
                    zscat2[i] = deltaR2 * round(tmp2/deltaR);            
                    }while((xscat2[i]*xscat2[i]+zscat2[i]*zscat2[i])>vmaxradscat2quad);
                }
            }//only sum2 if
        }//single and sum2 if
    else if(speckletype==spherical){
        zscat=(double*) malloc(Nscatterers*sizeof(double));
        double tmp1,tmp2, tmp3, vmaxradscatquad=vmaxradscat*vmaxradscat;
        for(int i=0;i<Nscatterers;i++){
            do{
                tmp1 = 2.0*vmaxradscat*frand()-vmaxradscat;
                tmp2 = 2.0*vmaxradscat*frand()-vmaxradscat;
                tmp3 = 2.0*vmaxradscat*frand()-vmaxradscat;
                //putting points in grid
                //Here round remplaces ANINT in fortran code
                xscat[i] = deltaR * round(tmp1/deltaR);
                yscat[i] = deltaR * round(tmp2/deltaR);
                zscat[i] = deltaR * round(tmp3/deltaR);                
                }while((xscat[i]*xscat[i]
                        +yscat[i]*yscat[i]
                        +zscat[i]*zscat[i])>vmaxradscatquad);
            }
        }//spherical if
    else if(speckletype==shell){
        zscat=(double*) malloc(Nscatterers*sizeof(double));
        double angle, urandom, urandom2, tmp1, tmp2, tmp3;
        for(int i=0;i<Nscatterers;i++){
            angle = 2.0*M_PI*frand()-vmaxradscat;
            urandom = 2.0*frand()-1.0;
            urandom2=urandom*urandom;

            tmp1 = sqrt(1.0 - urandom2)*cos(angle) * vmaxradscat;
            tmp2 = sqrt(1.0 - urandom2)*sin(angle) * vmaxradscat;
            tmp3 = urandom*vmaxradscat;
                
            //putting points in grid
            //Here round remplaces ANINT in fortran code
            xscat[i] = deltaR * round(tmp1/deltaR);
            yscat[i] = deltaR * round(tmp2/deltaR);
            zscat[i] = deltaR * round(tmp3/deltaR);
            }
        }
    else{
        fprintf(stderr,"Error: No definition for this speckletype call 1\n");
        abort();
        }

    //Create Speckle
    for(int nx3=0, x1=0, x2=0; nx3<npx; nx3++, x1+=npx2, x2+=npu2){
        if(speckletype==sum2){
            double temp=focal+dx*nx3;
            xcord[2]=temp;
            invlambdafoc2=1.0/(vlambda*focal2);
            }
        else if(speckletype==single){
            double temp=focal+dx*nx3;
            xcord[2]=temp;
            invlambdafoc2=1.0/(vlambda*temp);
            }
        else if((speckletype==spherical) || (speckletype==shell)) xcord[2]=dx*nx3;
        else{
            fprintf(stderr,"Error: No definition for this speckletype call 2\n");
            abort();            
            }
        double xscata, yscata, arg1, arg2;
        for(int nx2=0, y1=0, y2=0; nx2<npx; nx2++, y1+=npx, y2+=npu){
            
            xcord[1]=(speckletype==sum2?focal+dx*nx2:dx*nx2);
            
            for(int nx1=0, idx, idx2; nx1<npx; nx1++){
                idx=x1+y1+nx1;      //For all the arrays
                idx2=x2+y2+nx1;     //For VP that have extra size
                
                xcord[0]=dx*nx1;
                if(speckletype==single){
                    printf("%d %d %d\n", nx1, nx2,nx3);
                    double tmp;
                    for(int i=0;i<Nscatterers;i++){
                        xscata=xscat[i];
                        yscata=yscat[i];
                        
                        tmp=M_PI*(2.0*(xcord[0]*xscata+xcord[1]*yscata)   //arg
                                -(xscata*xscata+yscata*yscata)            //arg2
                                );    

                        arg1=tmp*invlambdafoc;
                        arg2=tmp*invlambdafoc2;
                        
                        ff3[idx]+=(
                                   alphascat[i]
                                   *cexp(betascat[i]*I)
                                   *(cexp(arg1*I)+cexp(arg2*I))
                                   );
                        }
                    }//single if
                else if(speckletype==sum2){
                    double tmp, xscata2, zscata2;
                    for(int i=0;i<Nscatterers;i++){
                        xscata=xscat[i];
                        yscata=yscat[i];
                        xscata2=xscat2[i];
                        zscata2=zscat2[i];      //yscata2 sustitutes zscata from fortran                        
                        arg1=M_PI*(2.0*(xcord[0]*xscata+xcord[1]*yscata)  //arg in Fortran
                                -(xscata*xscata+yscata*yscata)            //arg2 in Fortran
                                )*invlambdafoc;
                        
                        arg2=M_PI*(2.0*(xcord[0]*xscata2+xcord[1]*zscata2)   
                                -(xscata2*xscata2+zscata2*zscata2)      
                                )*invlambdafoc2;

                        ff3[idx]+=(
                                   alphascat[i]
                                   *cexp(betascat[i]*I)
                                   *cexp(arg1*I)
                                   );
                        
                        ff3_2[idx]+=(
                                     alphascat2[i]
                                     *cexp(betascat2[i]*I)
                                     *cexp(arg2*I)
                                     );
                        }
                    }//only sum2 if
                else if((speckletype==spherical)||(speckletype==shell)){
                    double zscata;
                    for(int i=0;i<Nscatterers;i++){
                        xscata=xscat[i];
                        yscata=yscat[i];
                        zscata=zscat[i];
                        arg1 =2.0*M_PI*invlambdafoc
                                 *(xcord[0]*xscata+xcord[1]*yscata+xcord[2]*zscata);
                        ff3[idx]+=(
                                   alphascat[i]
                                   *cexp(betascat[i]*I)
                                   *cexp(arg1*I)
                                   );
                        }
                    }//spherical and shell
                else{
                    fprintf(stderr,"Error: No definition for this speckletype\n");
                    abort();    
                    }
                if(speckletype!=sum2){
                    ffaux=ff3[idx];
                    //idx and idx2 are completly different
                    //and are calculated after the "for" calls
                    VP[idx2]=creal(ffaux*conj(ffaux));
                    aver+=VP[idx2];
                    }
                }//for nx1
            }//for nx2
        }//for nx3
    
    if(speckletype==sum2){
        double complex factor=(2.0*M_PI*dx/vlambda)*I;
        for(int nx3=0; nx3<npx; nx3++){
            for(int nx2=0; nx2<npx; nx2++){
                for(int nx1=0; nx1<npx; nx1++){
                    //remember nx2=npx*npx    
                    ffaux = ff3[npx2+nx2*npx+nx1]*cexp(rphase*I+factor*nx3)
                           +ff3_2[nx3*npx2+npx+nx1]*cexp(factor*nx2);
                    
                    int idx=nx3*npu2+nx2*npu+nx1;
                    VP[idx]=creal(ffaux*conj(ffaux));
                    aver+=VP[idx];
                    }//for nx1
                }//for nx2
            }//for nx1
        }//if sum2
    
    aver/=(npx2*npx);         //The average is calculated, remember npx2=npx*npx

    for(int nx3=0; nx3<npx; nx3++){
        double vaux;
        for(int nx2=0; nx2<npx; nx2++){
            for(int nx1=0, ivI; nx1<npx; nx1++){
                vaux=(VP[nx3*npu2+nx2*npu+nx1]/=aver);
                stdev+=(vaux*vaux);
                ivI=(int)(vaux/dvI);
                if((ivI<ngridvI)&&(ivI>=0)) pvI[ivI]+=adpp;
                }
            }
        }
    
    stdev=sqrt(stdev*adpp-1.0);
    printf("speckle st. dev.=%lf\n",stdev);
    
    if(speckletype==sum2){
        if(nrescale){
            printf("Shifting and rescaling field to have mean=1 and st.dev. = 1\n");
            double factor=1.0-1.0/stdev;
            //This moves also the unasigned values, but don't really matter
            for(int i=0; i<npu2*npu; i++){ 
                VP[i]=VP[i]/stdev+factor;
                }            
            }
        else{
            printf("NOT Shifting and rescaling field\n");
            }
        }
    

    const double factor=1.0-1.0/stdev;
    for(int i=0;i<ngridvI;i++){
        double vaux=(nrescale?(double(i)+0.5)*dvI/stdev+factor:((double)i+0.5)*dvI);
        fprintf(f18,"%d %lf %lf\n",i,vaux,pvI[i]/dvI);
        }

    //Fill end points (periodicity) in all the axis.
    //This is not the more efficient way, but more clear for to read
    //the iterator names were changed because are not fixed to any axis
    for(int i=0;i<npx;i++){
        for(int j=0;j<npx;j++){
            //In the same orther than in the Fortran code
            //The 0 is added only for readability
            VP[  i*npu2 +   j*npu + npx] = VP[ i*npu2 + j*npu + 0];
            VP[  i*npu2 + npx*npu +   j] = VP[ i*npu2 + 0*npu + j];
            VP[npx*npu2 +   i*npu +   j] = VP[ 0*npu2 + i*npu + j];
            }
        VP[  i*npu2 + npx*npu + npx] = VP[i*npu2];
        VP[npx*npu2 +   i*npu + npx] = VP[i*npu ];
        VP[npx*npu2 + npx*npu +   i] = VP[i     ];
        }
    VP[npx*npu2 + npx*npu + npx] = VP[0];
    
    //We can free memory here because this arrays are not used any more
    free(xscat);
    free(yscat);
    free(alphascat);
    free(betascat);
    if(speckletype==sum2){
        free(xscat2);
        free(zscat2);
        free(alphascat2);
        free(betascat2);
        free(ff3_2);
        }
    else if((speckletype==spherical)||(speckletype==shell)){
        free(zscat);
        }
    
    free(ff3);
    free(pvI);
    fclose(f18);
    printf("Speckle created\n");
    return 0;
    }


speckle::~speckle(){
    if(x) free(x);
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

void speckle::writespeckle(){
    const int npu=npxpu, npu2=npxpu*npxpu;
    
    printf("Writing Speckle\n");
    FILE *f14=fopen("speckle3D.dat","w");

    for(int i=0;i<npu;i++)
        for(int j=0;j<npu;j++)
            for(int k=0;k<npu;k++)
                fprintf(f14,"%lf\t %lf\t %lf\t %lf\n",xpos[k],xpos[j],xpos[i],VP[i*npu2+j*npu+k]);
    fclose(f14);
    printf("Speckle written\n");
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
    double *x_real=(double*) malloc(npx*npx2*sizeof(double));
    double complex *x_cmplx=(double complex*) malloc(padx*npx2*sizeof(double complex));
    
    if(Vk) free(Vk);
    Vk=(double complex *) malloc(npx*npx2*sizeof(double complex));

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
    status = DftiCreateDescriptor(&hand,DFTI_DOUBLE,DFTI_REAL,3,N);
    if(status!=0){
        printf("Error call DftiCreateDescriptor, status = %li\n",status);
        return status;
        }

    printf("Set out-of-place\n");
    status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if(status!=0){
        printf("Error 'placement' call DftiSetValue, status = %li\n",status);
        return status;
        }

    printf("Set CCE storage\n");
    status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if(status!=0){
        printf("Error 'event' call DftiSetValue, status = %li\n",status);
        return status;
        }

    status = DftiSetValue(hand, DFTI_FORWARD_SCALE, scaleforward);
    if(status!=0){
        printf("Error 'scale' call DftiSetValue, status = %li\n",status);
        return status;
        }

    status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rstrides);
    if(status!=0){
        printf("Error 'input stride' call DftiSetValue, status = %li\n",status);
        return status;
        }

    status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cstrides);
    if(status!=0){
        printf("Error 'input stride' call DftiSetValue, status = %li\n",status);
        return status;
        }

    printf("Commit DFTI descriptor\n");
    status = DftiCommitDescriptor(hand);
    if (status!=0){
        printf("Error call DftiCommitDescriptor, status = %li\n",status);
        return status;
        }

    printf("Compute forward transform\n");
    status = DftiComputeForward(hand, x_real, x_cmplx);
    if (status!=0){
        printf("Error call DftiComputeForward, status = %li\n",status);
        return status;
        }
                   
    printf("Reconfigure DFTI descriptor for backward transform\n");
    status = DftiSetValue(hand, DFTI_INPUT_STRIDES, cstrides);
    if(status!=0){
        printf("Error 'input stride' call DftiSetValue (back), status = %li\n",status);
        return status;
        }

    status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, rstrides);
    if(status!=0){
        printf("Error 'input stride' call DftiSetValue (back), status = %li\n",status);
        return status;
        }
    
    status = DftiSetValue(hand, DFTI_BACKWARD_SCALE, scalebackward);
    if(status!=0){
        printf("Error 'scale' call DftiSetValue (back), status = %li\n",status);
        return status;
        }

    status = DftiCommitDescriptor(hand);
    if (status!=0){
        printf("Error call DftiCommitDescriptor (back), status = %li\n",status);
        return status;
        }

    status = DftiComputeBackward(hand, x_cmplx, x_real);
    if (status!=0){
        printf("Error call DftiComputeBackward, status = %li\n",status);
        return status;
        }


    printf("Verify the result after a forward and backward Fourier Transforms\n");
    int count=0;
    for(int i=0;i<npx;i++){
        for(int j=0;j<npx;j++){
            for(int k=0;k<npx;k++){
                count+=(fabs(x_real[i*npx2+j*npx+k]-VP[i*npu2+j*npu+k])>0.01);
                }
            }
        }
    
    printf("Number of elements with an error greater than 0.01 %d\n",count);
    //This is completly modified way to do this because no mod or extra if are needed
    for(int i=0, ni=npx; i<npx; i++, ni--){
        for(int j=0, nj=npx; j<npx; j++, nj--){
            for(int k=0;k<padx;k++){ //same vales than the padded
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
    
    int* vkx=(int*)malloc(Ntot*sizeof(double)),
         vky=(int*)malloc(Ntot*sizeof(double)),
         vkz=(int*)malloc(Ntot*sizeof(double));

    double* diagTk=(double*)malloc(Ntot*sizeof(double));
    
    if(A) free(A);
    A=(double complex *) calloc(Ntot*Ntot,sizeof(double complex));
    
    printf("The dimension of the matrix is %d\n", Ntot);

    //----------Here we define the matrix A-----------------
    //vectors of the k-components and of the kinetic energy (in units of deltaT=hbar^2*unitk^2/2m)    
    for(int nkx=0,kx=-npmax/2,ntk=0 ;nkx<npmax; nkx++,kx++){         //ntk is declared here
        for(int nky=0,ky=-npmax/2 ;nky<npmax; nky++,ky++){
            for(int nkz=0,kz=-npmax/2; nkz<npmax; nkz++,kz++,ntk++){ //but only incremented here
                diagTk(ntk)=kx*kx+ky*ky+kz*kz;
                vkx(ntk)=kx;
                vky(ntk)=ky;
                vkz(ntk)=kz;
                }
            }
        }

    printf('Upper definition of A\n');
    //This loop is transverse to the cache access, but this is the way
    //it is made in the original code.
    //This can be made much more efficiently in the previous loop, but
    //maybe this way is usefull or readable. 
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
    }
