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

#include "specklemod.h"

int speckle::parse(){
    STARTDBG
    struct option longopts[] = {
    //first the variables available in the original code
    {"Nscatterers", required_argument, 0, 'N'},
    {"size", required_argument, 0, 's'},
    {"vDi", required_argument, 0, 'v'},
    {"vDi2", required_argument, 0, 'V'},
    {"focal", required_argument, 0, 'f'},
    {"focal2", required_argument, 0, 'F'},
    {"vlambda", required_argument, 0, 'l'},
    {"cofactor", required_argument, 0, 'c'},
    {"rphase", required_argument, 0, 'p'},
    {"npmax", required_argument, 0, 'P'},
    {"speckletype", required_argument, 0, 't'},
    {"nov", required_argument, 0, 'n'},
    {"intensity", required_argument, 0, 'i'},
    {"min", required_argument, 0, 'm'},      //min interestin value
    {"max", required_argument, 0, 'M'},      //max interestin value
    //now some extra variables for runtime
    {"begin", required_argument, 0, 'b'},    //starting point
    {"end", required_argument, 0, 'e'},      //end point
    {"solver", required_argument, 0, 'S'},   //solver type
    {"file", required_argument, 0, 'o'},     //output_prefix
    {"gpu", required_argument, 0, 'G'},      //ngpu
    {"restart",required_argument, 0, 'R'},   //continue previous calculations
    {"binsize",required_argument, 0, 'B'},   //histogram binsize
    {"dir",required_argument, 0, 'd'},       //directory for results
    {"save_interval",required_argument, 0, 'I'},       //directory for results
    
    {"gnuplot", no_argument, 0, 'g'},
    {"rescale", no_argument, 0, 'r'},
    {"vectors", no_argument, 0, 'a'},        
    {"help", no_argument, 0, 'h'},
        
    {0, 0, 0, 0}         //mandatory for end condition
        
    };

    size_t len=0;
    char *line, w1[50], w2[50];
    int c, option_index = 0, i=0, linenum=0;
    optind=1;    //this is a global variable reseted to restart reading the argumemts; 

    const char optstring[]="N:s:v:V:f:F:l:c:p:P:t:n:i:m:M:b:e:S:o:G:R:d:B:I:grah";
    
    // read all the options
    while((c=getopt_long(largc, largv, optstring,
                         longopts, &option_index))!=-1){
        if(c=='h'){
            if(rank==0){
                int c;
                FILE *file;
                file = fopen("help", "r");
                if (file) {
                    while ((c = getc(file)) != EOF)
                        putchar(c);
                    fclose(file);
                    }
                else{
                    fprintf(stderr,"no help file found\n");
                    }
                }
            exit(EXIT_SUCCESS);            
            }
        }
    if(optind>=largc){                       //error no Input file
        if(rank==0){
            fprintf(stderr,"Error: missing arguments, please check the command line\n");
            }
        exit(EXIT_FAILURE);
        }
    
    FILE* fp=fopen(largv[optind],"r");
    if(!fp){
        if(rank==0){
            fprintf(stderr,"Error input file \"%s\" couldn't be open\n",largv[optind]);
            }
        exit(EXIT_FAILURE);
        }
    
    //Parse input file
    while(getline(&line,&len,fp)>0){
        linenum++;
        if(line[0]!='#'){
            line=strtok(line,"#");
            if(sscanf(line,"%s %s %*s",w1,w2)>=1){
                i=0;
                while(true){
                    option curopt=longopts[i++];
                    if(curopt.name==0){
                        if(rank==0){
                            fprintf(stderr,"Warning: No valid option \"%s\" in line %d\n",w1,linenum);
                            }
                        exit(EXIT_FAILURE);
                        }
                    if(strcmp(w1,curopt.name)==0){
                        #ifdef DEBUG
                        fprintf(stderr,"Parsing: %s %s\n",w1,w2);
                        #endif //DEBUG
                        use_option(curopt.val,w2);
                        break;
                        }
                    }
                }
            }
        }
    fclose(fp);
    
    //parse command line arguments
    optind=1;
    while((c=getopt_long(largc, largv, optstring,
                         longopts, &option_index))!=-1){
        use_option(c,optarg);
        }
    ENDDBG
    return 0;
    }

void speckle::print(FILE* ou,char pre){
    fprintf(ou,"%c Nscatterers\t %d\n",pre,Nscatterers);
    fprintf(ou,"%c nrescale   \t %s\n",pre,nrescale?"true":"false");
    fprintf(ou,"%c gnuplot    \t %s\n",pre,usegnuplot?"true":"false");
    fprintf(ou,"%c size       \t %lg\n",pre,size);
    fprintf(ou,"%c vDi        \t %lg\n",pre,vDi);
    fprintf(ou,"%c vDi2       \t %lg\n",pre,vDi2);
    fprintf(ou,"%c focal      \t %lg\n",pre,focal);
    fprintf(ou,"%c focal2     \t %lg\n",pre,focal2);
    fprintf(ou,"%c vlambda    \t %g\n",pre,vlambda);    
    fprintf(ou,"%c cofactor   \t %lg\n",pre,cofactor);
    fprintf(ou,"%c rphase     \t %lg\n",pre,rphase);
    fprintf(ou,"%c npmax      \t %d\n",pre,npmax);
    fprintf(ou,"%c speckletype\t %s\n",pre,speckletype.c_str());
    fprintf(ou,"%c nov        \t %d\n",pre,nov);
    fprintf(ou,"%c vectors    \t %s\n",pre,vectors?"true":"false");
    fprintf(ou,"%c intensity  \t %lg%+lfI\n", pre, creal(Vintensity), cimag(Vintensity));

    fprintf(ou,"%c min        \t %lg\n",pre,min);
    fprintf(ou,"%c max        \t %lg\n",pre,max);
    fprintf(ou,"%c binsize    \t %lg\n",pre,binsize);

    fprintf(ou,"%c save_interv\t %d\n",pre,save_interval);
    
    fprintf(ou,"%c begin      \t %d\n",pre,start);
    fprintf(ou,"%c end        \t %d\n",pre,end);
                   
    fprintf(ou,"%c solver     \t %s\n",pre,solvername.c_str());
    fprintf(ou,"%c file_prefix\t %s\n",pre,fprefix.c_str());
    fprintf(ou,"%c save_dir   \t %s\n",pre,dirname.c_str());
    if(continuefile!=""){
        fprintf(ou,"%c continues  \t %s\n",pre,continuefile.c_str());
        }

    //Call the print from the parent
    base_calculator::print(ou,pre);
    }

void speckle::use_option(int opt,const char* thearg){
    switch (opt){
        case 'N': Nscatterers=atoi(thearg); break;
        case 's': size=atof(thearg); break;
        case 'v': vDi=atof(thearg); break;
        case 'V': vDi2=atof(thearg); break;
        case 'f': focal=atof(thearg); break;
        case 'F': focal2=atof(thearg); break;
        case 'l': vlambda=atof(thearg); break;
        case 'c': cofactor=atof(thearg); break;
        case 'p': rphase=atof(thearg); break;
        case 'P': npmax=atoi(thearg); break;
        case 't': speckletype=thearg; break;
        case 'n': nov=atoi(thearg); break;
        case 'r': nrescale=true; break;
        case 'a': vectors=true; break;
        case 'g': usegnuplot=true; break;
        case 'M': min=atof(thearg); break;
        case 'm': max=atof(thearg); break;
        case 'b': start=atoi(thearg); break;
        case 'e': end=atoi(thearg); break;
        case 'G': ngpu=atoi(thearg); break;
        case 'S': solvername=thearg; break;
        case 'o': fprefix=thearg; break;
        case 'R': continuefile=thearg; break;  
        case 'B': binsize=atof(thearg); break;
        case 'I': save_interval=atoi(thearg); break;
        case 'd': dirname=thearg; break;
            
        case 'i':
            double real, imag;
            sscanf(thearg,"%lf %lf[Ii]",&real,&imag);
            Vintensity=real+imag*I;
            break;
        case 'h': break;
        case ':':
        case '?':
            if(rank==0){
                fprintf(stderr,"Error: Wrong argument or option in command line\n");
                }
            exit(EXIT_FAILURE);
        }
    }

