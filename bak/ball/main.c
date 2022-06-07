#include <petsc.h>
#include <hdf5.h>
#include "plicit.h"

extern int Generator(int argc, char** argv);
extern int Explicit(int argc, char** argv);
extern int Implicit(int argc, char** argv);


int main(int argc, char** argv)
{   
    //~ generate the input data to HDF5 file
    if(strcmp(argv[1],"generator")==0)
        Generator(argc-1,argv+1);

    //~ iteration by explicit method
    else if(strcmp(argv[1],"explicit")==0)
        Explicit(argc-1,argv+1);

    //~ iteration by implicit method
    else if(strcmp(argv[1],"implicit")==0)
        Implicit(argc-1,argv+1);

    else
    {
        printf("You need to specify the function of this program.\n");
        return -1;
    }

    return 0;
}

