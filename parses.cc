#include "specklemod.h"

char Specklename[][10]={"spherical","sum2","single","shell"};

Speckletype atospeckle(const char* thename){
    int len=sizeof(Specklename)/sizeof(Specklename[0]);
    for(int i=0;i<len;i++){
        if(strcmp(thename,Specklename[i])==0){
            return Speckletype(i);
            }
        }
    fprintf(stderr,"Error invalid conversion from string to Speckletype\n");
    exit(EXIT_FAILURE);
    }

const char* speckletostr(Speckletype sp){
    int len=sizeof(Specklename)/sizeof(Specklename[0]);
    if((sp>=0)&&(sp<len)) return Specklename[sp];
    fprintf(stderr,"Error invalid conversion from Speckletype to string\n");
    exit(EXIT_FAILURE);
    }

void speckle::print_help(){
    printf("This is the help need to be updated\n");

    }

void speckle::print(FILE* ou){
    fprintf(ou,"Nscatterers\t %d\n",Nscatterers);
    fprintf(ou,"nrescale\t %s\n",nrescale?"true":"false");
    fprintf(ou,"size\t %lf\n",size);
    fprintf(ou,"vDi\t %lf\n",vDi);
    fprintf(ou,"focal\t %lf\n",focal);
    fprintf(ou,"vlambda\t %g\n",vlambda);
    fprintf(ou,"vDi2\t %lf\n",vDi2);
    fprintf(ou,"focal2\t %lf\n",focal2);
    fprintf(ou,"cofactor\t %lf\n",cofactor);
    fprintf(ou,"rphase\t %lf\n",rphase);
    fprintf(ou,"npmax\t %d\n",npmax);
    fprintf(ou,"npxpu\t %d\n",npxpu);
    fprintf(ou,"speckletype\t %s\n",speckletostr(speckletype));
    fprintf(ou,"nov\t %d\n",nov);
    }

void speckle::use_option(int opt,const char* thearg){
        switch (opt){
        case 'N':
            Nscatterers=atoi(thearg); break;
        case 's':
            size=atof(thearg); break;
        case 'v':
            vDi=atof(thearg); break;
        case 'V':
            vDi2=atof(thearg); break;
        case 'f':
            focal=atof(thearg); break;
        case 'F':
            focal2=atof(thearg); break;
        case 'l':
            vlambda=atof(thearg); break;
        case 'c':
            cofactor=atof(thearg); break;
        case 'p':
            rphase=atof(thearg); break;
        case 'm':
            npmax=atoi(thearg); break;
        case 't':
            speckletype=atospeckle(thearg); break;
        case 'n':
            nov=atoi(thearg); break;
        case 'r':
            nrescale=true; break;
        case 'h':
            print_help();
            exit(EXIT_SUCCESS);
        case ':':
        case '?':
            fprintf(stderr,"Error: Wrong argument or option in command line\n");
            exit(EXIT_FAILURE);
            }    
    }

void speckle::parser(){

    if (largc==1){
        printf("\tUsage: executable [-opt1 val1] [-opt2 val2]... input_file\n");
        printf("\ttype \"executable -h\" for the full list of valid arguments\n\n");
        exit(EXIT_FAILURE);
        }

    //this is a global variable reseted to restart reading the argumemts; 
    FILE* fp;
    size_t len=0;
    char *line, w1[50], w2[50];
    int c, option_index = 0, i=0, linenum=0;
    optind=1;
    
    //This looks terrible, I know, but it is the saver way
    while((c=getopt_long(largc, largv, "N:s:v:V:f:F:l:c:p:m:s:n:rh",
                         longopts, &option_index))!=-1){
        if(c=='h') use_option(c,optarg);
        };
    if(optind>=largc){                       //error no Input file
        fprintf(stderr,"Error input file, missing argument\n");
        exit(EXIT_FAILURE);
        }
    if((fp=fopen(largv[optind],"r"))==NULL){
        fprintf(stderr,"Error input file \"%s\" couldn't be open\n",largv[optind]);
        exit(EXIT_FAILURE);
        }
    //Parse input file
    while(getline(&line,&len,fp)>0){
        linenum++;
        if(line[0]!='#'){
            line=strtok(line,"#");
            if(sscanf(line,"%s %s %*s",w1,w2)>=1){
                int tmp_deb=0;
                i=0;
                while(true){
                    option curopt=longopts[i++];
                    if(curopt.name==0){
                        fprintf(stderr,"Warning: No valid option \"%s\" in line %d\n",
                                w1,linenum);
                        break;
                        }
                    if(strcmp(w1,curopt.name)==0){
                        use_option(curopt.val,w2);
                        break;
                        }
                    if(++tmp_deb==20) break;
                    }
                }
            }
        }
    
    //parse command line arguments
    optind=1;
    while((c=getopt_long(largc, largv, "N:s:v:V:f:F:l:c:p:m:s:n:rh",
                        longopts, &option_index))!=-1){
        use_option(c,optarg);
        }
    }
