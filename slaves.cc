#include "slaves.h"

slaves::slaves(int argc, char** argv):
    speckle(argc, argv),
    cont(0),
    total_processed(0),    
    seeds(NULL),
    seeds_processed(NULL),
    logfile(NULL)
{
    STARTDBG
    #ifdef UMPI
    MPI_Init (&argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &wsize ); 
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    #endif
    if(rank==0){
        seeds=new int[wsize]();
        seeds_processed=new int[wsize]();
        if(!seeds){
            fprintf(stderr,"Error: seeds arrays allocation failed\n");
            printme();
            exit(EXIT_FAILURE);
            }
        pthread_mutex_init(&mutex1, NULL);
        pthread_mutex_init(&mutex2, NULL);
        pthread_mutex_init(&mutex3, NULL);
        
        //start logfile
        string filename= save_dir+fprefix+"logfile_"+speckletype+"_"+timestr+".log";
        
        logfile= fopen(filename.c_str(), "w");
        log_printf("Logfile for %s version of speckle code\n",VERSION);
        }
    ENDDBG
    }


void slaves::run(){
    STARTDBG
    int value_get, info;
    if(rank==0){
        value_get=getseed(0,0);
        #ifdef UMPI                 // this is interesting
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        // The function "Manager" sends the initial seeds                    
        pthread_create(&thread, &attr,
                       slaves::master_static, (void *)this);
        log_printf("Master thread created\n");
        }
    else{
        MPI_Recv(&value_get, 1, MPI_INT, 0, seedtag, MPI_COMM_WORLD, &status2);
        #endif
        }
    while(value_get!=-1){
        //Now start making calculations and check that no errors ocurred
        info=calculate(value_get);
        int rec=(info==-1?-1:indices);
        if(rank==0) {
            value_get=getseed(0,rec);
            #ifdef UMPI
            }
        else{
            MPI_Sendrecv(&rec, 1, MPI_INT, 0, resultag,
                         &value_get, 1,MPI_INT, 0, seedtag,
                         MPI_COMM_WORLD, &status2);
            #endif // UMPI
            }
        // if error calculating at this point rec==value_get==-1
        if(rec>0){
            if(rank==0) process(rec,values);
            #ifdef UMPI
            else{
                MPI_Send(values, rec, MPI_DOUBLE, 0, valuestag, MPI_COMM_WORLD);
                }
            #endif // UMPI
            }
        }
    #ifdef UMPI
    if(rank==0){
        pthread_join(thread, NULL);
        log_printf("pthread_join out\n");        
        }
    #endif // UMPI
    ENDDBG
    }

#ifdef UMPI
int slaves::master(){
    STARTDBG
    int rec=0;               // received value
    double *tmp_buffer;      // temporal buffer
    int source;              // sender process            

    for(int i=1;i<wsize;i++){
        int tmp=getseed(i,rec);
        MPI_Send(&tmp, 1, MPI_INT, i, seedtag, MPI_COMM_WORLD);
        }

    log_printf("Remote processes %d\n",cont);
    
    while(cont>0){
        MPI_Recv(&rec, 1,MPI_INT, MPI_ANY_SOURCE,
                 resultag, MPI_COMM_WORLD, &status);
        source=status.MPI_SOURCE;
        int seed=getseed(source,rec);
        // Send the new seed if no error, or -1 as confirmation
        MPI_Send(&seed, 1, MPI_INT, source, seedtag, MPI_COMM_WORLD);
        if(rec>0){
            tmp_buffer=(double*) malloc(rec*sizeof(double));
            MPI_Recv(tmp_buffer, rec, MPI_DOUBLE, source,
                     valuestag, MPI_COMM_WORLD, &status);
            process(rec,tmp_buffer);
            free(tmp_buffer);
            }
        else{
            printf("Thread_Warning: Got %d from process %d\n",rec,source);
            }
        log_printf("Processed %d results from process %d\n",rec,source);
        }
    log_printf("Finish thread\n");
    ENDDBG
    return 0;
    }

#endif // UMPI
