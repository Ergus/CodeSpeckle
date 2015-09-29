#include "specklemod.h"
#include "parallel.h"

int main(int argc, char **argv){
    printf("Test speckle\n");
    slaves a(argc,argv);
    printf("Running \n");
    a.run();

    return 0;
    }
