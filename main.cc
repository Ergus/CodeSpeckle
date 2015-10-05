#include "specklemod.h"
#include "slaves.h"

int main(int argc, char **argv){
    printf("Running \n");
    slaves a(argc,argv);
    a.run();
    printf("End Succesfull!!!\n");
    return 0;
    }
