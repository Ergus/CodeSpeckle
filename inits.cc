#include "specklemod.h"

int speckle::init_spherical(int idseed){
    const int npx=npmax, npx2=npx*npx, npu=npxpu,
        npu2=npxpu*npxpu, ngridvI =1000;
    
    double vmaxradscat =1.0/(2*vDi ),
        deltaR =vlambda*focal /size,
        vmaxradscatquad=vmaxradscat*vmaxradscat,
        aver=0.0,
        invlambdafoc=1.0/(vlambda*focal),
        stdev=0.0,
        vImax=10.0,
        dvI=vImax/ngridvI,
        adpp=1.0/(npx2*npx),
        xcord[3];

    double complex ffaux;             //this is used as a temporal variable, but in different loops    
    
    if(strcmp(Specklename[speckletype],"spherical")!=0){
        fprintf(stderr,"Error: call %s but Specklename: %s\n",
                __PRETTY_FUNCTION__, Specklename[speckletype]);
        return 1;
        }
    printf("Speckletype spherical\n");
    
    double *xscat=(double*) malloc(Nscatterers*sizeof(double)),
        *yscat=(double*) malloc(Nscatterers*sizeof(double)),
        *alphascat=(double*) malloc(Nscatterers*sizeof(double)),
        *betascat=(double*) malloc(Nscatterers*sizeof(double)),
        *zscat=(double*) malloc(Nscatterers*sizeof(double)),
        *pvI=(double*) calloc(ngridvI,sizeof(double));
    
    //This will have dimension 3
    double complex *ff3=(double complex*) calloc(npmax*npmax*npmax,sizeof(double complex));

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

    //sample positions of random scatterers
    for(int i=0;i<Nscatterers;i++){
        do{
            double tmp1 = 2.0*vmaxradscat*frand()-vmaxradscat;
            double tmp2 = 2.0*vmaxradscat*frand()-vmaxradscat;
            double tmp3 = 2.0*vmaxradscat*frand()-vmaxradscat;
            //putting points in grid
            //Here round remplaces ANINT in fortran code
            xscat[i] = deltaR * round(tmp1/deltaR);
            yscat[i] = deltaR * round(tmp2/deltaR);
            zscat[i] = deltaR * round(tmp3/deltaR);                
            }while((xscat[i]*xscat[i]
                    +yscat[i]*yscat[i]
                    +zscat[i]*zscat[i])>vmaxradscatquad);
        }

    //Create Speckle
    for(int nx3=0, x1=0, x2=0; nx3<npx; nx3++, x1+=npx2, x2+=npu2){
        xcord[2]=dx*nx3;
        for(int nx2=0, y1=0, y2=0; nx2<npx; nx2++, y1+=npx, y2+=npu){
            xcord[1]=dx*nx2;
            for(int nx1=0, idx, idx2; nx1<npx; nx1++){
                idx=x1+y1+nx1;      //For all the arrays
                idx2=x2+y2+nx1;     //For VP that have extra size
                xcord[0]=dx*nx1;
                for(int i=0;i<Nscatterers;i++){
                    double xscata=xscat[i];
                    double yscata=yscat[i];
                    double zscata=zscat[i];
                    double arg1 =2.0*M_PI*invlambdafoc
                        *(xcord[0]*xscata+xcord[1]*yscata+xcord[2]*zscata);
                    ff3[idx]+=(alphascat[i]
                               *cexp(betascat[i]*I)
                               *cexp(arg1*I) );
                    }
                ffaux=ff3[idx];
                //idx and idx2 are completly different
                //and are calculated after the "for" calls
                VP[idx2]=creal(ffaux*conj(ffaux));
                aver+=VP[idx2];                    
                }//for nx1
            }//for nx2
        }//for nx3

    //The average is calculated, remember npx2=npx*npx
    aver/=(npx2*npx);

    for(int nx3=0; nx3<npx; nx3++){
        for(int nx2=0; nx2<npx; nx2++){
            for(int nx1=0, ivI; nx1<npx; nx1++){
                double vaux=(VP[nx3*npu2+nx2*npu+nx1]/=aver);
                stdev+=(vaux*vaux);
                ivI=(int)(vaux/dvI);
                if((ivI<ngridvI)&&(ivI>=0)) pvI[ivI]+=adpp;
                }
            }
        }

    stdev=sqrt(stdev*adpp-1.0);
    printf("speckle st. dev.=%lf\n",stdev);
    
    const double factor=1.0-1.0/stdev;
    for(int i=0;i<ngridvI;i++){
        double vaux=(nrescale?(double(i)+0.5)*dvI/stdev+factor:((double)i+0.5)*dvI);
        fprintf(f18,"%d %lf %lf\n",i,vaux,pvI[i]/dvI);
        }

    //Fill end points (periodicity) in all the axis.
    //This is not the most efficient way, but more clear to read
    //the iterator names were changed because they are not fixed to any axis
    for(int i=0;i<npx;i++){
        for(int j=0;j<npx;j++){
            //In the same orther than in the Fortran code
            //The 0 is added only for readability
            VP[  i*npu2 +  j*npu +npx] = VP[i*npu2 +j*npu +0];
            VP[  i*npu2 +npx*npu +  j] = VP[i*npu2 +0*npu +j];
            VP[npx*npu2 +  i*npu +  j] = VP[0*npu2 +i*npu +j];
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
    free(zscat);
    
    free(ff3);
    free(pvI);
    fclose(f18);
    printf("Speckle created\n");
    return 0;
    }

int speckle::init_sum2(int idseed){
    //For to be used internally, not needed but for to make
    //the code looks like the Fortran One
    const int npx=npmax, npx2=npx*npx, npu=npxpu, npu2=npxpu*npxpu;
    
    const int ngridvI =1000;
    double *pvI=(double*) calloc(ngridvI,sizeof(double));

    double vmaxradscat =1.0/(2*vDi ),
        vmaxradscat2=1.0/(2*vDi2),
        deltaR =vlambda*focal /size,
        deltaR2=vlambda*focal2/size,        
        vmaxradscatquad=vmaxradscat*vmaxradscat,
        vmaxradscat2quad=vmaxradscat2*vmaxradscat2,
        invlambdafoc=1.0/(vlambda*focal),
        invlambdafoc2,               //this is only references in plain and sum2
        aver=0.0,        
        stdev=0.0,
        vImax=10.0,
        dvI=vImax/ngridvI,
        adpp=1.0/(npx2*npx),
        xcord[3];
    
    //this is used as a temp, but in different loops
    double complex ffaux;

    if(strcmp(Specklename[speckletype],"sum2")!=0){
        fprintf(stderr,"Error: call %s but Specklename: %s\n",
                __PRETTY_FUNCTION__, Specklename[speckletype]);
        return 1;
        }
    printf("Speckletype sum2\n");    
    
    double *xscat=(double*) malloc(Nscatterers*sizeof(double)),
        *yscat=(double*) malloc(Nscatterers*sizeof(double)),
        *alphascat=(double*) malloc(Nscatterers*sizeof(double)),
        *betascat=(double*) malloc(Nscatterers*sizeof(double)),
        *xscat2=(double*) malloc(Nscatterers*sizeof(double)),
        *zscat2=(double*) malloc(Nscatterers*sizeof(double)),
        *alphascat2=(double*) malloc(Nscatterers*sizeof(double)),
        *betascat2=(double*) malloc(Nscatterers*sizeof(double));

    //This will have dimension 3
    double complex *ff3=(double complex*) calloc(npmax*npmax*npmax,sizeof(double complex)),
        *ff3_2=(double complex*) calloc(npmax*npmax*npmax,sizeof(double complex));        

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
        alphascat2[i] = 2.0*frand();
        betascat2[i] = 2.0*M_PI*frand()-M_PI;       
        }
    
    //sample positions of random scatterers
    for(int i=0;i<Nscatterers;i++){
        do{
            double tmp1 = 2.0*vmaxradscat*frand()-vmaxradscat;
            double tmp2 = 2.0*vmaxradscat*frand()-vmaxradscat;
            //putting points in grid
            //Here round remplaces ANINT in fortran code
            xscat[i] = deltaR * round(tmp1/deltaR);
            yscat[i] = deltaR * round(tmp2/deltaR);
            }while((xscat[i]*xscat[i]+yscat[i]*yscat[i])>vmaxradscatquad);
        do{
            double tmp1 = 2.0*vmaxradscat2*frand()-vmaxradscat2;
            double tmp2 = 2.0*vmaxradscat2*frand()-vmaxradscat2;
            //putting points in grid
            //Here round remplaces ANINT in fortran code
            xscat2[i] = deltaR2 * round(tmp1/deltaR);
            zscat2[i] = deltaR2 * round(tmp2/deltaR);         
            }while((xscat2[i]*xscat2[i]+zscat2[i]*zscat2[i])>vmaxradscat2quad);
        }

    //Create Speckle
    for(int nx3=0, x1=0, x2=0; nx3<npx; nx3++, x1+=npx2, x2+=npu2){
        double temp=focal+dx*nx3;
        xcord[2]=temp;
        invlambdafoc2=1.0/(vlambda*focal2);
        for(int nx2=0, y1=0, y2=0; nx2<npx; nx2++, y1+=npx, y2+=npu){
            xcord[1]=focal+dx*nx2;
            for(int nx1=0, idx, idx2; nx1<npx; nx1++){
                idx=x1+y1+nx1;      //For all the arrays
                idx2=x2+y2+nx1;     //For VP that have extra size                
                xcord[0]=dx*nx1;
                for(int i=0;i<Nscatterers;i++){
                    double xscata=xscat[i];
                    double yscata=yscat[i];
                    double xscata2=xscat2[i];
                    double zscata2=zscat2[i];        //yscata2 sustitutes zscata from fortran
                    double arg1=M_PI*(2.0*(xcord[0]*xscata+xcord[1]*yscata)  //arg in Fortran
                               -(xscata*xscata+yscata*yscata)                //arg2 in Fortran
                               )*invlambdafoc;
                        
                    double arg2=M_PI*(2.0*(xcord[0]*xscata2+xcord[1]*zscata2)   
                               -(xscata2*xscata2+zscata2*zscata2)      
                               )*invlambdafoc2;

                    ff3[idx]+=(alphascat[i]
                               *cexp(betascat[i]*I)
                               *cexp(arg1*I));
                        
                    ff3_2[idx]+=(alphascat2[i]
                                 *cexp(betascat2[i]*I)
                                 *cexp(arg2*I));
                    }
                }//for nx1
            }//for nx2
        }//for nx3
    
    const double complex aux=(2.0*M_PI*dx/vlambda)*I;
    for(int nx3=0; nx3<npx; nx3++){
        for(int nx2=0; nx2<npx; nx2++){
            for(int nx1=0; nx1<npx; nx1++){
                //remember nx2=npx*npx    
                ffaux = ff3[npx2+nx2*npx+nx1]*cexp(rphase*I+aux*nx3)
                    +ff3_2[nx3*npx2+npx+nx1]*cexp(aux*nx2);
                    
                int idx=nx3*npu2+nx2*npu+nx1;
                VP[idx]=creal(ffaux*conj(ffaux));
                aver+=VP[idx];
                }//for nx1
            }//for nx2
        }//for nx1
    
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
    const double factor=1.0-1.0/stdev;
    
    printf("speckle st. dev.=%lf\n",stdev);
    
    if(nrescale){
        printf("Shifting and rescaling field to have mean=1 and st.dev. = 1\n");
        //This moves also the unasigned values, but don't really matter
        for(int i=0; i<npu2*npu; i++){ 
            VP[i]=VP[i]/stdev+factor;
            }       
        }
    else{
        printf("NOT Shifting and rescaling field\n");
        }
    
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
    free(xscat2);
    free(zscat2);
    free(alphascat2);
    free(betascat2);
    free(ff3_2);

    free(ff3);
    free(pvI);
    fclose(f18);
    printf("Speckle created\n");
    return 0;
    }



int speckle::init_single(int idseed){
    //For to be used internally, not needed but for to make
    //the code looks like the Fortran One
    const int npx=npmax, npx2=npx*npx, npu=npxpu,
        npu2=npxpu*npxpu, ngridvI =1000;

    double vmaxradscat =1.0/(2*vDi ),
        vmaxradscat2=1.0/(2*vDi2),
        deltaR =vlambda*focal /size,
        vmaxradscatquad=vmaxradscat*vmaxradscat,
        invlambdafoc=1.0/(vlambda*focal),
        aver=0.0,
        stdev=0.0,
        vImax=10.0,
        dvI=vImax/ngridvI,
        adpp=1.0/(npx2*npx),
        xcord[3];

    double complex ffaux;             //this is used as a temp, but in different loops    
    
    if(strcmp(Specklename[speckletype],"single")!=0){
        fprintf(stderr,"Error: call %s but Specklename: %s\n",
                __PRETTY_FUNCTION__, Specklename[speckletype]);
        return 1;
        }
    printf("Speckletype single\n");
    
    double *xscat=(double*) malloc(Nscatterers*sizeof(double)),
        *yscat=(double*) malloc(Nscatterers*sizeof(double)),
        *alphascat=(double*) malloc(Nscatterers*sizeof(double)),
        *betascat=(double*) malloc(Nscatterers*sizeof(double)),
        *pvI=(double*) calloc(ngridvI,sizeof(double));

    //This will have dimension 3
    double complex *ff3=(double complex*) calloc(npmax*npmax*npmax,sizeof(double complex));

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

    //sample positions of random scatterers
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

    //Create Speckle
    for(int nx3=0, x1=0, x2=0; nx3<npx; nx3++, x1+=npx2, x2+=npu2){
        double temp=focal+dx*nx3;
        xcord[2]=temp;
        double invlambdafoc2=1.0/(vlambda*temp);
        for(int nx2=0, y1=0, y2=0; nx2<npx; nx2++, y1+=npx, y2+=npu){
            
            xcord[1]=(speckletype==sum2?focal+dx*nx2:dx*nx2);
            for(int nx1=0, idx, idx2; nx1<npx; nx1++){
                idx=x1+y1+nx1;      //For all the arrays
                idx2=x2+y2+nx1;     //For VP that have extra size
                
                xcord[0]=dx*nx1;
                for(int i=0;i<Nscatterers;i++){
                    double xscata=xscat[i];
                    double yscata=yscat[i];
                    
                    double tmp=M_PI*(2.0*(xcord[0]*xscata+xcord[1]*yscata) //arg
                              -(xscata*xscata+yscata*yscata)               //arg2
                              );

                    double arg1=tmp*invlambdafoc;
                    double arg2=tmp*invlambdafoc2;
                        
                    ff3[idx]+=(alphascat[i]
                               *cexp(betascat[i]*I)
                               *(cexp(arg1*I)+cexp(arg2*I))
                               );
                    }
                
                ffaux=ff3[idx];
                //idx and idx2 are completly different
                //and are calculated after the "for" calls
                VP[idx2]=creal(ffaux*conj(ffaux));
                aver+=VP[idx2];

                }//for nx1
            }//for nx2
        }//for nx3
    
    aver/=(npx2*npx);         //The average is calculated, remember npx2=npx*npx

    for(int nx3=0; nx3<npx; nx3++){
        for(int nx2=0; nx2<npx; nx2++){
            for(int nx1=0, ivI; nx1<npx; nx1++){
                double vaux=(VP[nx3*npu2+nx2*npu+nx1]/=aver);
                stdev+=(vaux*vaux);
                ivI=(int)(vaux/dvI);
                if((ivI<ngridvI)&&(ivI>=0)) pvI[ivI]+=adpp;
                }
            }
        }
    
    stdev=sqrt(stdev*adpp-1.0);
    printf("speckle st. dev.=%lf\n",stdev);
    
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
    free(ff3);
    free(pvI);
    fclose(f18);
    printf("Speckle created\n");
    return 0;
    }

int speckle::init_shell(int idseed){
    //For to be used internally, not needed but for to make
    //the code looks like the Fortran One
    const int npx=npmax, npx2=npx*npx, npu=npxpu,
        npu2=npxpu*npxpu, ngridvI =1000;
    
    double vmaxradscat =1.0/(2*vDi ),
        vmaxradscat2=1.0/(2*vDi2),
        deltaR =vlambda*focal /size,
        deltaR2=vlambda*focal2/size,        
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
    
    double complex ffaux;             //this is used as a temp, but in different loops    

    if(strcmp(Specklename[speckletype],"shell")!=0){
        fprintf(stderr,"Error: call %s but Specklename: %s\n",
                __PRETTY_FUNCTION__, Specklename[speckletype]);
        return 1;
        }
    printf("Speckletype shell\n");    
    
    double *pvI=(double*) calloc(ngridvI,sizeof(double)),
        *xscat=(double*) malloc(Nscatterers*sizeof(double)),
        *yscat=(double*) malloc(Nscatterers*sizeof(double)),
        *alphascat=(double*) malloc(Nscatterers*sizeof(double)),
        *betascat=(double*) malloc(Nscatterers*sizeof(double)),
        *zscat=(double*) malloc(Nscatterers*sizeof(double));
    
    //This will have dimension 3
    double complex *ff3=(double complex*) calloc(npmax*npmax*npmax,sizeof(double complex));

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

    for(int i=0;i<Nscatterers;i++){
        double angle = 2.0*M_PI*frand()-vmaxradscat;
        double urandom = 2.0*frand()-1.0;
        double urandom2=urandom*urandom;
        
        double tmp1 = sqrt(1.0-urandom2)*cos(angle)*vmaxradscat;
        double tmp2 = sqrt(1.0-urandom2)*sin(angle)*vmaxradscat;
        double tmp3 = urandom*vmaxradscat;
                
        //putting points in grid
        //Here round remplaces ANINT in fortran code
        xscat[i] = deltaR * round(tmp1/deltaR);
        yscat[i] = deltaR * round(tmp2/deltaR);
        zscat[i] = deltaR * round(tmp3/deltaR);
        }

    //Create Speckle
    for(int nx3=0, x1=0, x2=0; nx3<npx; nx3++, x1+=npx2, x2+=npu2){
        xcord[2]=dx*nx3;
        
        for(int nx2=0, y1=0, y2=0; nx2<npx; nx2++, y1+=npx, y2+=npu){    
            xcord[1]=dx*nx2;
            
            for(int nx1=0, idx, idx2; nx1<npx; nx1++){
                idx=x1+y1+nx1;      //For all the arrays
                idx2=x2+y2+nx1;     //For VP that have extra size
                xcord[0]=dx*nx1;

                for(int i=0;i<Nscatterers;i++){
                    double xscata=xscat[i];
                    double yscata=yscat[i];
                    double zscata=zscat[i];
                    double arg1 =2.0*M_PI*invlambdafoc
                        *(xcord[0]*xscata+xcord[1]
                          *yscata+xcord[2]*zscata);
                    ff3[idx]+=(alphascat[i]
                               *cexp(betascat[i]*I)
                               *cexp(arg1*I) );
                    }

                ffaux=ff3[idx];
                //idx and idx2 are completly different
                //and are calculated after the "for" calls
                VP[idx2]=creal(ffaux*conj(ffaux));
                aver+=VP[idx2];
                }//for nx1
            }//for nx2
        }//for nx3
    
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
            VP[  i*npu2 +  j*npu +npx] =VP[i*npu2 +j*npu +0];
            VP[  i*npu2 +npx*npu +  j] =VP[i*npu2 +0*npu +j];
            VP[npx*npu2 +  i*npu +  j] =VP[0*npu2 +i*npu +j];
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
    free(zscat);
    free(ff3);
    free(pvI);
    fclose(f18);
    printf("Speckle created\n");
    return 0;
    }
