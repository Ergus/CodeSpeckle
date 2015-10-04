
#ifndef BASE_CALCULATOR_H
#define BASE_CALCULATOR_H

#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <time.h>

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

using namespace std;

//Template elements are <derived class, type of produced data, number of arrays>

/********************************************//**
 * \brief Usage for the base class defined for allow master-slave using templates.
 *
 * This class is a base class template needed to implement easy paralelizable 
 * using master-slave thread bases system. In general the class defines a virtual
 * int calculate(), an array of pointers to functions for procesing the results
 * and some extra variables for mpi management
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
        base_calculator():
            start(1),                      
            end(1),                        
            rank(0),
            wsize(1),
            local_rank(0),
            local_size(1),
            ngpu(0),      // you should set this value in the input file
            ncpu(sysconf(_SC_NPROCESSORS_ONLN)-1),
            hostname(getenv("HOSTNAME")){
            //Initialize time variables just one time at very begining
            time (&rawtime);
            timeinfo = localtime (&rawtime);
            strftime(timestr,80,"%a_%F_%X",timeinfo);
            }

        /// Get the arrays with the dimensions
        /** \return Array of dimensions for every result. */
        int get_indices(){return indices;}

        /// Get the arrays with the results.
        /** \return The results array number i result[i].
            \note i must be lower than the dimension of the problem. */
        double* get_values(){return values;};


        /// Calculate method, this is very important.
        /** \param [in] seed interger seed for the random number generator.
            \return The return value should be 0, else an error ocurred
           
            This method should make all the most important calculations.
            The results should be stored (or just direct the pointer) in the 
            "results" array. And the size (number of results) stored in the 
            #indices. */
        virtual int calculate(int seed)=0;

        /// Process function for results.
        /** \param [in] nvalues Interger number (size) of results in this group.
            \param [in] ovalues Array with results.
            \return The return value should be 0, else an error ocurred
           
            This is a base class method, but it really calls methods defined in the 
            derived class. should be set. */
        virtual int process(int nvalues, double* ovalues)=0;

    protected:
        /** \name Run control variables
            Run system check needs this variables to start and stop
            in this implementation a serial for is used, but this can be reimplemented 
            in the future. */                
        ///\{
        int start,                ///< First seed
            end;                  ///< Last seed
        ///\}

        /** \name MPI variables
            This are needed mainly for error control system. In not MPI compilation, 
            then the default values provided in this constructor are the right for 
            the rest of the program execution. */
        ///\{
        int rank,                 ///< MPI rank in world
            wsize,                ///< MPI World size
            local_rank,           ///< Rank in local node
            local_size;           ///< Number of processes running in local node
        int ngpu,                 ///< GPU number
            ncpu;                 ///< CPU number
        ///\}
        
        string hostname;          ///< Name of the host (current node)
        int indices;              ///< Number of elements
        double *values;           ///< Array with the results

        /** \name Some time variables */
        ///\{
        time_t rawtime;           ///< time_t object for start time
        struct tm * timeinfo;     ///< tm struct for start time
        char timestr[80];         ///< char* with start time %a_%F_%X
        ///\}
    };

#endif // BASE_CALCULATOR_H
