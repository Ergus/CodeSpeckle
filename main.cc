#include "specklemod.h"
#include "parser.h"

int main(int argc, char **argv){
    printf("Test speckle\n");

    parser par(argc,argv);
    
    printf("Creating \n");
    speckle a(&par);
    printf("Running \n");
    a.init(1);
    printf("Transforming \n");
    a.ftspeckle();
    printf("Defining \n");
    a.defineA();

    return 0;
    }
