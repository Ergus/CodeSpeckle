#include "specklemod.h"

int main(int argc, char **argv){
    printf("Test speckle\n");

    
    printf("Creating \n");
    speckle a(argc,argv);
    a.init(100);
    a.ftspeckle();
    a.defineA();

    return 0;
    }
