#include "specklemod.h"

int main(){
    printf("Test speckle\n");

    
    printf("Creating Single\n");
    speckle a(10,sum2);
    a.init_sum2(100);
    a.ftspeckle();
    a.defineA();
    
    return 0;
    }
