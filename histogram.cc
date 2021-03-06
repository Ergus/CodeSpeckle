/******************************************************
 * This file is part of the Codespeckle++ distribution Copyright (c) 2015 Jimmy
 * Aguilar Mena.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************/

#include "histogram.h"

histogram::histogram(speckle *outter):
    deltaE(outter->binsize),
    save_interval(outter->save_interval),
    nrealiz(0),
    dos(NULL),sumr(NULL),meanr(NULL),nos(NULL),countr(NULL),
    gnuplot(NULL),caller(outter),filename("")
{
    STARTDBG;
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

    // this is to define the new output only, this can change easily
    filename=caller->dirname+"/"+caller->filename+"_values.dat";
    
    f10=fopen(filename.c_str(),"w");       // values.dat    

    // Check all the pointers are fine
    if (!(nos && countr && dos && sumr && meanr && f10)){
            fprintf(stderr,"Error: allocating memory or opening file for histogram\n");
            printme();
            exit(EXIT_FAILURE);
            }
    
    // Continue previous calculations if indicated
    if (caller->continuefile!=""){
        load_previous(caller->continuefile.c_str());
        }
    if (caller->usegnuplot){
        if(caller->wsize==1){
            gnuplot=popen("gnuplot -persistent", "w");
            if (!gnuplot){
                fprintf(stderr,"Gnuplot option active, but pipe couldn't be open\n");
                }
            }
        else{
            fprintf(stderr,"Use gnuplot option only in serial code\n");
            }
        }
    
    ENDDBG;
    }

histogram::~histogram(){
    caller->log_printf("Saving values at the end\n");
    write_values();
    free(nos);
    free(countr);            
    free(dos);
    free(sumr);
    free(meanr);
                    
    fclose(f10);
    if(gnuplot) pclose(gnuplot);
    }

int histogram::load_previous(const char *filename){
    STARTDBG;
        FILE *fp=fopen(filename,"r");
    if (!fp){
        fprintf(stderr,"Impossible import file %s to continue calculations\n",
                filename);
        return -1;
        }
    char ignore[1024];
    fgets(ignore, sizeof(ignore), fp);  //ignore first 2 lines
    fgets(ignore, sizeof(ignore), fp);
    //start parsing
    for(int i=0;i<ndeltaE;i++){
        int matched=fscanf(fp,"%*d %*f %lf %lf %lf %d",
                           &dos[i], &meanr[i], &sumr[i], &countr[i]);

        // an error exit the program because some
        // values maybe were modified.
        if (matched<4){  
            fprintf(stderr,"Error importing file %s line %d\n",
                    filename,i+2);
            exit(EXIT_FAILURE);
            }
        }
    fclose(fp);
    ENDDBG;
    return(0);
    }

int histogram::process(const int np,const double *values){
    STARTDBG;
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

    if((nrealiz%save_interval)==0){
        caller->log_printf("Saving values nrealiz=%d\n",nrealiz);
        write_values();
        }

    ENDDBG;
    return 0;
    }

int histogram::write_values(){
    STARTDBG;
    rewind(f10);
    fprintf(f10,"#Emax= %lf, deltaE= %lf, nrealiz=%d\n",
            EMax,deltaE,nrealiz);
    fprintf(f10,"#n E(n) dos(E) meanr(E) sumr(E) countr(E)\n");
                    
    double E=0;
    for(int i=0;i<ndeltaE;i++){
        E=(i+1)*deltaE;
        meanr[i]=(countr[i]==0 ? 0.0 : sumr[i]/countr[i]);
        fprintf(f10,"%d %lf %lf %lf %lf %d\n",
                i, E, dos[i], meanr[i], sumr[i], countr[i]);
        }
    fflush(f10);
        
    if (gnuplot){
        printf("plotting process %d\n",caller->rank);
        fprintf(gnuplot,"plot \"%s\" u 2:5\n",filename.c_str());
        fflush(gnuplot);
        }
    ENDDBG;
    return 0;
    }
