
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>

#ifdef UMPI
#include <mpi.h>
namespace internal{
    inline MPI_Datatype mpi_type_id (const int *){return MPI_INT;}
    inline MPI_Datatype mpi_type_id (const long int *){return MPI_LONG;}
    inline MPI_Datatype mpi_type_id (const unsigned int *){return MPI_UNSIGNED;}
    inline MPI_Datatype mpi_type_id (const unsigned long int *){return MPI_UNSIGNED_LONG;}
    inline MPI_Datatype mpi_type_id (const unsigned long long int *){return MPI_UNSIGNED_LONG_LONG;}
    inline MPI_Datatype mpi_type_id (const float *){return MPI_FLOAT;}
    inline MPI_Datatype mpi_type_id (const double *){return MPI_DOUBLE;}
    inline MPI_Datatype mpi_type_id (const long double *){return MPI_LONG_DOUBLE;}
    }

#define resultag 0     // Tag for results
#define seedtag 1      // Tag for seeds
#define valuestag 1    //

#endif

#include "parser.h"

// to reuse this code compile with and without -DUMPI.

using namespace std;

template<typename CALCULATOR, typename T,int NGROUPS>
class slaves{
    public:
        slaves(parser *thepar):
            rank(0), // My id in the worst case always 0 
            size(1), // The world size, again at least always 1
            begin(thepar->start),
            end(thepar->end),
            it(thepar->start),
            cont(0),
            seeds(NULL),
            mutex2(NULL)
            {
                #ifdef UMPI
                MPI_Init (&thepar->argc, &thepar->argv);
                MPI_Comm_size( MPI_COMM_WORLD, &size ); 
                MPI_Comm_rank( MPI_COMM_WORLD, &rank );
                MPI_T=internal::mpi_type_id(&T);
                #endif
                if(rank==0){
                    seeds=new int[size];
                    mutex2=new int[size];
                    }
                calc=new CALCULATOR(thepar);
                }

        
        ~slaves(){
            delete calc;
            if(seeds) delete seeds;  //only deallocate in rank 0
            if(mutex2) delete mutex2;
            #ifdef UMPI
            MPI_Finalize();
            #endif
            }

        int run(){
            // here the names are the same but the work compleetly different            
            int value_get, info;
            int rec[NGROUPS]={};     // received value
            T *tmp_buffer;                        
            if(rank==0){
                value_get=get_seed(0,rec);
                #ifdef UMPI  //this is interesting
                pthread_attr_init(&attr);
                pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                // The function "Manager" sends the initial seeds                    
                pthread_create(&thread, &attr,
                               slaves::master_static, (void *)&this); 
                }
            else{
                MPI_Recv(&value_get, 1, MPI_INT, 0, seedtag, MPI_COMM_WORLD, &status2);
                #endif
                }
            while(value_get!=-1){
                info=calc->run(value_get);
                rec=calc->get_rets();
                if(rank==0) {
                    value_get=getseed(0,rec);
                    #ifdef UMPI
                    }
                else{
                    MPI_Sendrecv(rec, NGROUPS, MPI_INT, 0, resultag,
                                 &value_get,  1, MPI_INT, 0, seedtag,
                                 MPI_COMM_WORLD, &status2);                    
                    #endif // UMPI
                    }
                
                for(int i=0;i<NGROUPS;i++){
                    if(rec[i]!=0){
                        if(rank==0) process(i,rec[i],calc->results[i]);
                        #ifdef UMPI
                        else{
                            MPI_Send(calc->results[i], rec[i], MPI_T, 0,
                                     valuestag, MPI_COMM_WORLD);
                            #endif // UMPI
                            }
                        }
                    }
                }
            }

        void process(int index,int nvalues, T* array){
            pthread_mutex_lock(&mutex2[index]);        // mutex for threads save
            calc->process[index](nvalues,array);
            pthread_mutex_unlock(&mutex2[index]);      // mutex for threads save
            }

        #ifdef UMPI //       
        // this will load the master as a daemon, because a
        // pthread cannot receive a normal function
        static void* master_static(void *in){
            slaves *in2=(slaves *) in;
            in2->master();
            pthread_exit((void *)NULL);
            }

        // this is the real master
        int master(){
            int rec[NGROUPS]={};     // received value
            T *tmp_buffer;            
            int source;              // sender process            
            MPI_Status status;

            for(int i=1;i<size;i++){
                MPI_Send(getseed(i,rec), 1, MPI_INT, i, seedtag, MPI_COMM_WORLD);
                }

            while(cont>0){
                MPI_Recv(rec, NGROUPS,MPI_INT, MPI_ANY_SOURCE,
                         resultag, MPI_COMM_WORLD, &status);
                
                source=status.MPI_SOURCE;
                seed=input->getseed(source,rec);
                if(rec[0]>-1){
                    // Send the new seed if no error
                    MPI_Send(&seed, 1, MPI_INT, source, seedtag, MPI_COMM_WORLD);
                    if(rec[0]>0){
                        for(int i=0;i<NGROUPS;i++){
                            if(rec[i]!=0){
                                tmp_buffer=(double*) malloc(rec[i]*sizeof(T));
                                MPI_Recv(tmp_buffer, rec[i], MPI_T, source,
                                         valuestag, MPI_COMM_WORLD, &status);
                                
                                process(i,rec[i],tmp_buffer);
                                free(tmp_buffer);
                                }
                            }
                        }
                    }
                else{
                    printf("Thread_Warning: Got %d from process %d\n",rec,source);
                    }
                }
            }
        #endif // UMPI
        
        int* getseed(int proc, int *returned){
            pthread_mutex_lock(&mutex1);      // mutex for threads save
            
            if((proc!=0)&&(seeds[proc]==0)) cont++;
            if((it<=end)&&(returned[0]!=-1)){
                seeds[proc]=it++;
                }
            else{
                if(returned[0]==-1){
                    fprintf(stderr,"Error: Process %d failed and sent %d, last seed send to it was: %d\n", proc,returned,seeds[proc]);
                    }
                seeds[proc]=-1;
                if (proc!=0) cont--;
                }
            // unlock mutex, the thread is protected by mpi, and if it is called
            // localy don't cares
            pthread_mutex_unlock(&mutex1);    
            return &seeds[proc];
            }
        
    private:
        int sizes[NGROUPS];     //will get the size of every individual array
        T **arrays;             //will get the pointers arrays
        int rank, size;         //mpi constants
        CALCULATOR *calc;
        
        #ifdef UMPI             
        pthread_t thread;        // thread will be loaded in the paralell version
        pthread_attr_t attr;     // Attribute for thread
        MPI_Datatype MPI_T;      // Datatype templatized in constructor
        #endif

        //static members 
        int it, cont;       //it is iterator
        const int begin, end; //start end values (including bounds)
        int *seeds;
        pthread_mutex_t mutex1,    // Mutex for getseed
            *mutex2;    // array of mutex for process results

        
    };
