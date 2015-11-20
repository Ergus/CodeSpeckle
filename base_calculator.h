
#ifndef BASE_CALCULATOR_H
#define BASE_CALCULATOR_H

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <pthread.h>     // obvious
#include <signal.h>      // for threads mutex
#include <getopt.h>      // get enviroment variables
#include <stdarg.h>      // for the personalized print functions
#include <string>
#include <string.h>   //strtok and strcmp in parser

#include <sys/sysinfo.h> // get_nprocq
#include <sys/stat.h>    // mkdir

#ifdef UMPI
#include <mpi.h>
#define resultag 0     // Tag for results
#define seedtag 1      // Tag for seeds
#define valuestag 1    // Tag for values
#define VERSION "parallel"
#else
#define VERSION "serial"
#endif

//This is a debbuger macro for fft functions that needs to
//return 0 each speckle will print its own error message before interrupt

#ifndef dbg
#define dbg(x) {                                                           \
        if(x!=0){                                                          \
            fprintf(stderr,"Error: %s returned %d\n",#x, x);               \
            printme();                                                     \
            return(-1);                                                    \
            }                                                              \
        }
#endif  //dbg

#ifndef dbg_mem
#define dbg_mem(x) {                                                       \
        if(x==NULL){                                                       \
            fprintf(stderr,"Error: %s is %p, after allocation\n",#x, x);     \
            printme();                                                     \
            return(-1);                                                    \
            }                                                              \
        }
#endif  //dbg_mem

#ifndef printme
#ifdef MPI_VERSION
#define printme() {                                                  \
        fprintf(stderr,"%s in %s:%d (process %s in %s)\n",           \
                __PRETTY_FUNCTION__, __FILE__, __LINE__,             \
                getenv("OMPI_COMM_WORLD_RANK"), getenv("HOSTNAME")); \
        }
#else  //MPI_VERSION
#define printme() {                                        \
        fprintf(stderr,"%s in %s:%d\n",                    \
                __PRETTY_FUNCTION__, __FILE__, __LINE__);  \
        }
#endif // MPI_VERSION
#endif // printme

#ifdef DEBUG
 #define STARTDBG fprintf(stderr,"Calling %s: %s:%d\n",                   \
                         __PRETTY_FUNCTION__,__FILE__,__LINE__);
 #define ENDDBG fprintf(stderr,"Ended %s: %s:%d\n",                       \
                         __PRETTY_FUNCTION__,__FILE__,__LINE__);
#else
 #define STARTDBG
 #define ENDDBG
#endif

using namespace std;

//Template elements are <derived class, type of produced data, number of arrays>

/********************************************//**
 * \brief Usage for the base class defined for allow master-slave.
 *
 * This class is a base class template needed to implement easy paralelizable 
 * using master-slave thread bases system. In general the class defines a virtual
 * methods for calculating and procesing the results.
 * The pure virtual methods #getseed_serial, #process_serial and #calculate are 
 * mandatory to be implemented in any derived class.
 *
 * \author Jimmy Aguilar Mena
 * \version 0.1
 ***********************************************/
class base_calculator{
    public:
        /// Default constructor for base class.
        /** No parameters needed. Any modification should be inserted latter
            in the derived class constructor, not here please.
            You don't need to declare anything here, this is a template. */
        base_calculator(int argc, char **argv);

        /// Destructor for the base calculator.
        ~base_calculator();
        
        /// Get the arrays with the dimensions
        /** \return Array of dimensions for every result. */
        inline int get_indices(){return indices;}

        /// Get the arrays with the results.
        /** \return The results array number i result[i].
            \note i must be lower than the dimension of the problem. */
        inline double* get_values(){return values;};

        /// Run routine, is the manager for all the system
        /** It calls the #process, #getseed and #calculate routines
            that must be implementedin any derived class. */
        int run();

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
            fprintf(logfile,"[ %s %6ld ] ",buffer,lt-rawtime);
            vfprintf(logfile, format, args);
            va_end(args);
            fflush(logfile);
            pthread_mutex_unlock(&mutex3);
            }
        

        /// Print the values that will be used in this class.
        /** \param [in,out] ou to print, default standard output
            \param [in] pre char to put before every line, default '#'
            The function needs to be defined in the derived class and call 
            #print_parent.*/
        virtual void print(FILE* ou=stdout,char pre='#');
        
    protected:
        /// Initialization for the parallel framework (or serial).
        /** Here are set some initializations that can not be defined in the
            constructor because some variables are defined in the derived class 
            so this should be called at the end or after the derived 
            class constructor. Now is inserted in the begin of #run so no explicit
            call*/
        int initialize();

        /** \name Virtual methods
            The next 3 methods are pure virtual and need to be defined in the derived
            class implementations. all are mandatory. */
        ///\{
        /// Serial version for the getseed routine
        /** This function needs to be defined in any derived class. 
            Because the paralell version needs it it is an overloaded function, 
            but the other one is called only one calls this */
        virtual int getseed_serial()=0;

        /// Process function for results.
        /** \param [in] nvalues Interger number (size) of results in this group.
            \param [in] ovalues Array with results.
            \return The return value should be 0, else an error ocurred
           
            This is needs to be defined in the derived class with the function
            that will process the results. */
        virtual int process_serial(int nvalues, double* ovalues)=0;
        
        /// Calculate method, this is very important.
        /** \param [in] seed interger seed for the random number generator.
            \return The return value should be 0, else an error ocurred
           
            This method should make all the most important calculations.
            The results should be stored (or just direct the pointer) in the 
           "results" array. And the size (number of results) stored in the 
            #indices. */
        virtual int calculate(int seed)=0;        
        ///\}
        
        /** \name MPI variables
            This are needed mainly for error control system. In not MPI compilation, 
            then the default values provided in this constructor are the right for 
            the rest of the program execution. */
        ///\{
        int rank,                     ///< MPI rank in world
            wsize,                    ///< MPI World size
            local_rank,               ///< Rank in local node
            local_size;               ///< Number of processes running in local node
        int ngpu,                     ///< GPU number
            ncpu,                     ///< CPU number/process
            start_cpu;                ///< Indef for the first process affinity in this process
        ///\}

        /** \name Usefull strings
            This is the name (general pattern) for the output file and some 
            other usefull strings. */
        ///\{
        string hostname;              ///< Name of the host (current node).
        string filename;              ///< General name for the files generated.
        string dirname;               ///< General name for the directory name.
        string timestr;               ///< String with start time %a_%F_%X
        ///\}

        /** \name Result variables 
            the results for the calculations should be set in this variables*/
        ///\{
        int indices;                  ///< Number of elements.
        double *values;               ///< Array with the results.
        ///\}
        
        /** \name Some time variables */
        ///\{
        time_t rawtime;               ///< time_t object for start time
        struct tm * timeinfo;         ///< tm struct for start time
        ///\}

        int largc;                    ///< Standard input counter
        char **largv;                 ///< Standard input array
        
    private:
        #ifdef UMPI                    
        pthread_t thread;             ///< Thread will be loaded in the paralell version
        pthread_attr_t attr;          ///< Attribute for thread
        MPI_Status status, status2;   ///< tha status for the mpi processes
        cpu_set_t cpus;
        #endif

        //static members 
        int cont,                     ///< Cont is the number of remote processes.
            total_processed;          ///< Number of total processed systems
        int *seeds,                   ///< Seeds in every remote
            *seeds_processed;         ///< Count of succesfull processed systems
        pthread_mutex_t mutex1,       ///< Mutex for getseed
            mutex2,                   ///< Mutex for process function
            mutex3;                   ///< Mutex for log_printf

        FILE *logfile;                ///< Logfile pointer

        /// Routine for process results. (thread save)
        /** This is a thread save call for the inheritated #process serial. 
            This is called only by the run() function*/
        int process(int nvalues, double* array);

        /// Complex routine to get the seeds
        /** The original code is implemented serialy, but this one 
            is thread save and with MPI support. */
        int getseed(int proc, int returned);

        #ifdef UMPI
        /// master static helper class.
        /** \param in object of type slave pointer. 
            This is a helper function, because pthread don't accepts member functions
            This will load the master as a daemon, because a
            pthread cannot receive a normal function. */
        static void* master_static(void *in){
            STARTDBG
            base_calculator *in2=(base_calculator *) in;
            in2->master();
            pthread_exit((void *)NULL);
            ENDDBG
            }

        /// This is the real master function.
        /** To be called by the pthreads objects it needs the helper function 
            #master_static. Because pthreads don't accepts not static member functions.
            but non static members don't have access to members.*/
        int master();
        
        #endif // UMPI                
        
    };

#endif // BASE_CALCULATOR_H
