#include "specklemod.h"

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
    nov(4){
        VP=(double*) malloc(npxpu*npxpu*npxpu*sizeof(double));
        xpos=(double*) malloc(npxpu*sizeof(double));
        x=(double*) malloc(npmax*sizeof(double));

        dx=size/npmax;
        dxi=double(npmax)/size;
        }

void speckle::init(int idseed){
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
        //fprintf(f18,"%d %lf %lf\n",i,vaux,pvI[i]/dvI);
        fprintf("%d %lf %lf\n",i,vaux,pvI[i]/dvI);
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
    }


speckle::~speckle(){
    if(x) free(x);
    if(xpos) free(xpos);    
    if(VP) free(VP);
    }
