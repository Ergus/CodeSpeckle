#include "parser.h"

parser::parser(int argc,char** argv,
               int orank, int oworld_size):
    //Constructor using standard input.
    Nscatterers(100),
    npmax(0),                 //This values are checked for errors
    nov(4),    
    size(0.11),
    vDi(0.05635245901639344262),
    vDi2(0.05635245901639344262),    
    focal(40.0),
    focal2(30.0),
    vlambda(800.0e-6),
    cofactor(1.0),    
    rphase(0.0),
    nrescale(false),    	  //shift and rescale sum2 speckle        
    //run needed
    start(0), fprefix("output"),
    //MPI needed
    rank(orank), world_size(oworld_size),
    local_rank(0), local_size(1),
    largc(argc), largv(argv){
    
    hostname=getenv("HOSTNAME");

    if (largc==1){
        if(rank==0){
            printf("\tUsage: executable [--opt1 val1] [--opt2 val2]... input_file\n");
            printf("\ttype \"executable -h\" for the full list of valid arguments\n\n");
            }
        exit(EXIT_FAILURE);
        }
    
    if(rank>0){
        local_rank = atoi( getenv("OMPI_COMM_WORLD_LOCAL_RANK"));
        local_size = atoi( getenv("OMPI_COMM_WORLD_LOCAL_SIZE"));
        }
    
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
        {"npmax", required_argument, 0, 'm'},
        {"speckletype", required_argument, 0, 't'},
        {"nov", required_argument, 0, 'n'},
        {"rescale", no_argument, 0, 'r'},
        
        //now some extra variables for runtime
        {"start", required_argument, 0, '1'},    //starting point
        {"end", required_argument, 0, '2'},      //end point
        {"solver", required_argument, 0, 'S'},   //solver type
        {"file", required_argument, 0, 'P'},     //output_prefix

        //please improve this option as needed
        {"help", no_argument, 0, 'h'},
        
        //mandatory for end condition        
        {0, 0, 0, 0}
        };

    size_t len=0;
    char *line, w1[50], w2[50];
    int c, option_index = 0, i=0, linenum=0;
    optind=1;    //this is a global variable reseted to restart reading the argumemts; 
    
    //This looks terrible, I know, but it is the saver way
    while((c=getopt_long(largc, largv, "N:s:v:V:f:F:l:c:p:m:s:n:rh",
                         longopts, &option_index))!=-1){
        if(c=='h') use_option(c,optarg);
        };
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
                int tmp_deb=0;
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
                        use_option(curopt.val,w2);
                        break;
                        }
                    if(++tmp_deb==20) break;
                    }
                }
            }
        }
    fclose(fp);
    
    //parse command line arguments
    optind=1;
    while((c=getopt_long(largc, largv, "N:s:v:V:f:F:l:c:p:m:s:n:rh",
                         longopts, &option_index))!=-1){
        use_option(c,optarg);
        }    

    }

void parser::print_help(){
    if(rank==0){
        printf("This is the help needs to be updated\n");
        }
    }

void parser::print(FILE* ou,char pre){
    fprintf(ou,"%c Nscatterers\t %d\n",pre,Nscatterers);
    fprintf(ou,"%c nrescale\t %s\n",pre,nrescale?"true":"false");
    fprintf(ou,"%c size\t %lf\n",pre,size);
    fprintf(ou,"%c vDi\t %lf\n",pre,vDi);
    fprintf(ou,"%c vDi2\t %lf\n",pre,vDi2);
    fprintf(ou,"%c focal\t %lf\n",pre,focal);
    fprintf(ou,"%c focal2\t %lf\n",pre,focal2);
    fprintf(ou,"%c vlambda\t %g\n",pre,vlambda);    
    fprintf(ou,"%c cofactor\t %lf\n",pre,cofactor);
    fprintf(ou,"%c rphase\t %lf\n",pre,rphase);
    fprintf(ou,"%c npmax\t %d\n",pre,npmax);
    fprintf(ou,"%c speckletype\t %s\n",pre,speckletype.c_str());
    fprintf(ou,"%c nov\t %d\n",pre,nov);


    fprintf(ou,"%c start\t %d\n",pre,start);
    fprintf(ou,"%c end\t %d\n",pre,end);
                   
    fprintf(ou,"%c solver\t %s\n",pre,solver.c_str());
    fprintf(ou,"%c file prefix\t %s\n",pre,fprefix.c_str());
    
    }

void parser::use_option(int opt,const char* thearg){
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
            speckletype=thearg; break;
        case 'n':
            nov=atoi(thearg); break;
        case 'r':
            nrescale=true; break;
        case '1':
            start=atoi(thearg); break;
        case '2':
            end=atoi(thearg); break;
        case 'S':
            solver=thearg; break;
        case 'P':
            fprefix=thearg; break;
        case 'h':
            print_help();
            exit(EXIT_SUCCESS);
        case ':':
        case '?':
            if(rank==0){
                fprintf(stderr,"Error: Wrong argument or option in command line\n");
                }
            exit(EXIT_FAILURE);
        }
    }
