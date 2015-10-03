#ifndef SLAVE_H
#define SLAVE_H

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include <string>
#include <time.h>

#ifdef UMPI
#include <mpi.h>
#define resultag 0     // Tag for results
#define seedtag 1      // Tag for seeds
#define valuestag 1    // Tag for values
#define VERSION "parallel"
#else
#define VERSION "serial"
#endif


#include "specklemod.h"

// to reuse this code compile with and without -DUMPI.
// De default makefile will do it for you

using namespace std;

/********************************************//**
* \brief Usage for the slave class
* 
* This is the main class that can be used easily for a thread based
* master-slave system.
* This class can work as manager for any class derived from #base_calculator
* that have defined a #getseed(), #calculate() and #process() routines
*
* \author Jimmy Aguilar Mena
* \version 0.1
***********************************************/
class slaves: public speckle{
    public:

        /// Constructor for class slave
        /** \param [in] argc standar input counter
            \param [in] argv standar input list of arguments. */
        slaves(int argc, char** argv);

        ~slaves(){
            STARTDBG
            pthread_mutex_destroy(&mutex1);
            pthread_mutex_destroy(&mutex2);
            pthread_mutex_destroy(&mutex3);
            if(seeds) delete[] seeds;  //only deallocate in rank 0
            if(seeds_processed) delete[] seeds_processed;
            if (logfile){
                log_printf("Slaves class destructed in process: %d\n",rank);
                fclose(logfile);
                }
            #ifdef UMPI
            MPI_Finalize();
            #endif
            ENDDBG
            }

        /// Run routine, is the manager for all the system
        /** It calls the inheritated routines (maybe reimplemented)
         * #getseed and #calculate. */
        void run();

        /// Routine for process results. (thread save)
        /** This is a thread save call for the inheritated #process routine
            if the original routine is thread save it is not needed to reimplement it. */
        int process(int nvalues, double* array){
            STARTDBG
            pthread_mutex_lock(&mutex2);        // mutex for threads save
            speckle::process(nvalues,array);
            pthread_mutex_unlock(&mutex2);      // mutex for threads save
            ENDDBG
            return 0;
            }

#ifdef UMPI
        /// master static helper class.
        /** \param in object of type slave pointer. 
            This is a helper function, because pthread don't accepts member functions
            This will load the master as a daemon, because a
            pthread cannot receive a normal function. */
        static void* master_static(void *in){
            STARTDBG
            slaves *in2=(slaves *) in;
            in2->master();
            pthread_exit((void *)NULL);
            ENDDBG
            }

        /// This is the real master function.
        /** To be called by the pthreads objects it needs the helper function 
            #master_static. Because pthreads don't accepts not static member functions.
            but non static members don't have access to members.
        */
        int master();
        
#endif // UMPI

        /// Complex routine to get the seeds
        /** The original code is implemented serialy, but this one 
            is thread save and with MPI support. */
        int getseed(int proc, int returned){
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
                tmp=speckle::getseed();          // get the new seed
                seeds_processed[proc]++;         // count the processed seeds 
                }
            if ((proc!=0) && (tmp==-1)) cont--;
            log_printf("Process %d returned %d results for seed %d "
                       "local processed %d total processed %d sending now %d\n",
                       proc,returned,seeds[proc],seeds_processed[proc],
                       total_processed, tmp);
            seeds[proc]=tmp;
            pthread_mutex_unlock(&mutex1);
            ENDDBG
            return seeds[proc];
            }
        
    protected:
        /// printf function that saves information in the logfile
        /** The losfile name is automatically defined in the constructor of 
            the slaves in the serial and paralell version for the code.
            Using this function the format for the logfile is fixed and saves 
            always the time of the event.
            Only the master process will save information in the log file.
            The arguments are the same than printf. */
        void log_printf(const char * format, ... ){
            pthread_mutex_lock(&mutex3);
            time_t lt;
            char buffer[80];
            time (&lt);
            struct tm * timeinfo = localtime(&lt);
            strftime(buffer,80,"%a %d/%m/%y %X",timeinfo);
    
            va_list args;
            va_start(args, format);
            fprintf(logfile,"[%s] ",buffer);
            vfprintf(logfile, format, args);
            va_end(args);
            fflush(logfile);
            pthread_mutex_unlock(&mutex3);
            }
        
    private:
#ifdef UMPI                    
        pthread_t thread;             // thread will be loaded in the paralell version
        pthread_attr_t attr;          // Attribute for thread
        MPI_Status status, status2;
#endif

        //static members 
        int cont,                     // cont is the number of remote processes running.
            total_processed;
        int *seeds, *seeds_processed;
        pthread_mutex_t mutex1,       // Mutex for getseed
            mutex2,                   // array of mutex for process results
            mutex3;

        FILE *logfile;
        time_t rawtime;
        struct tm * timeinfo;
        
    };

#endif //SLAVE_H
