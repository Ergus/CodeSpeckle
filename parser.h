#ifndef parser_h
#define parser_h

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <getopt.h>
#include <algorithm>

using namespace std;

class parser{
    public:
        parser(int argc,char** argv,
               int orank=0, int oworld_size=1);
        
        //Speckle variables
        int Nscatterers, npmax, nov;
        double size, vDi, vDi2, focal, focal2,
            vlambda, cofactor, rphase;
        bool nrescale;
        string speckletype;

        //Run important variables
        int start, end;
        string solver, fprefix;

        //mpi needed variables
        int rank, world_size;
        int local_rank, local_size;
        string hostname;
        
        //now some functions
        void print(FILE* ou=stdout,char pre=' ');

        void print_help();
        void use_option(int opt, const char* thearg);
    private:
        int largc;
        char **largv;
        
    };

#endif
