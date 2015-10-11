#include "base_calculator.h"

base_calculator::base_calculator(int argc, char **argv):
    rank(0),
    wsize(1),
    local_rank(0),
    local_size(1),
    ngpu(0),      // you should set this value in the input file
    ncpu(get_nprocs()),
    hostname(getenv("HOSTNAME")),
    filename(""),
    dirname("outputs"),
    largc(argc), largv(argv),
    cont(0),
    total_processed(0),
    seeds(NULL),
    seeds_processed(NULL),
    logfile(NULL){

    STARTDBG;
    
    // Start logfile
    #ifdef UMPI
    MPI_Init (&largc, &largv);
    MPI_Comm_size( MPI_COMM_WORLD, &wsize ); 
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    local_rank=atoi(getenv("OMPI_COMM_WORLD_LOCAL_RANK"));
    local_size=atoi(getenv("OMPI_COMM_WORLD_LOCAL_SIZE"));
    int tmpncpu=ncpu/local_size;
    ncpu=tmpncpu+(local_rank<(ncpu-tmpncpu*local_size));
    if(rank==0) ncpu--; //This is to let a free processor for the thread
    #endif    
    
    // Initialize time variables just one time at very begining
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    char t[80];
    strftime(t,80,"%a_%F_%X",timeinfo);
    timestr=t;
    
    ENDDBG;
    }

base_calculator::~base_calculator(){
    STARTDBG;
    pthread_mutex_destroy(&mutex1);
    pthread_mutex_destroy(&mutex2);
    pthread_mutex_destroy(&mutex3);
    if(seeds) delete[] seeds;  //only deallocate in rank 0
    if(seeds_processed) delete[] seeds_processed;
    if (logfile){
        time_t time_end;
        time(&time_end);
        log_printf("Elapsed time with np %d %d\n",wsize,time_end-rawtime);
        log_printf("Slaves class destructed in process: %d\n",rank);
        fclose(logfile);
        }
    #ifdef UMPI
    MPI_Finalize();
    #endif
    ENDDBG;
    }

int base_calculator::initialize(){
    STARTDBG;
    if(rank==0){
        // Prepare memory for master
        seeds=new int[wsize](); dbg_mem(seeds);
        seeds_processed=new int[wsize](); dbg_mem(seeds_processed);
        
        pthread_mutex_init(&mutex1, NULL);
        pthread_mutex_init(&mutex2, NULL);
        pthread_mutex_init(&mutex3, NULL);

        int created=mkdir(dirname.c_str(),0777);
        if (created==0) printf("Creating dir %s\n",dirname.c_str());

        //Define logfile
        string logfilename= dirname+"/"+filename+".log";
        
        logfile= fopen(logfilename.c_str(), "w");
        log_printf("Logfile for %s version of speckle code\n",VERSION);
        
        print(logfile);
        }
    ENDDBG;
    return 0;
    }

void base_calculator::run(){
    STARTDBG;
    
    int value_get=-1, info=initialize();

    if (rank==0){        
        value_get=getseed(0,0);
        #ifdef UMPI                 // this is interesting
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        // The function "Manager" sends the initial seeds                    
        pthread_create(&thread, &attr,
                       base_calculator::master_static, (void *)this);
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
    ENDDBG;
    }


int base_calculator::getseed(int proc, int returned){
    STARTDBG
        pthread_mutex_lock(&mutex1);      // mutex for threads save
    int tmp;
    // first check if is initial seed
    if(seeds[proc]==0) cont+=(proc!=0);
    else if(returned!=-1) total_processed++;
    // then check for errors
    if(returned==-1){
        fprintf(stderr,"Error: Process %d failed, sent %d, last seed was: %d\n",
                proc,returned,seeds[proc]);
        tmp=-1;
        }
    else{
        tmp=getseed_serial();            // get the new seed
        seeds_processed[proc]++;         // count the processed seeds 
        }
    if ((proc!=0) && (tmp==-1)) cont--;
    log_printf("Process %d returned %d results for seed %d "
               "local processed %d total processed %d sending now %d\n",
               proc,returned,seeds[proc],seeds_processed[proc],
               total_processed, tmp);
    seeds[proc]=tmp;
    pthread_mutex_unlock(&mutex1);
    ENDDBG;
    return seeds[proc];
    }

int base_calculator::process(int nvalues, double* array){
    STARTDBG;
    pthread_mutex_lock(&mutex2);        // mutex for threads save
    process_serial(nvalues,array);
    pthread_mutex_unlock(&mutex2);      // mutex for threads save
    ENDDBG;
    return 0;
    }

void base_calculator::print(FILE* ou,char pre){
    fprintf(ou,"%c start_time \t %s\n",pre,timestr.c_str());
    fprintf(ou,"%c rank       \t %d\n",pre,rank);
    fprintf(ou,"%c wsize      \t %d\n",pre,wsize);
    fprintf(ou,"%c ngpu       \t %d\n",pre,ngpu);
    fprintf(ou,"%c ncpu       \t %d\n",pre,ncpu);
    fprintf(ou,"%c filename   \t %s\n",pre,filename.c_str());
    }

#ifdef UMPI
int base_calculator::master(){
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
