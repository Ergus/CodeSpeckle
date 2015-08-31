#include "specklemod.h"

int main(){
    printf("Test speckle\n");

    
    printf("Creating \n");
    speckle a(10,spherical);
    a.init(100);
    a.ftspeckle();
    a.defineA();

    printf("Creating sum2\n");
    speckle b(10,sum2);
    b.init(100);
    b.ftspeckle();
    b.defineA();

    printf("Creating single\n");
    speckle c(10,single);
    c.init(100);
    c.ftspeckle();
    c.defineA();

    printf("Creating shell\n");
    speckle d(10,shell);
    d.init(100);
    d.ftspeckle();
    d.defineA();    
    
    return 0;
    }
