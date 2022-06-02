#include <petsc.h>
void H5W(int argc, char** argv);
int main(int argc, char** argv)
{   
    //  // 打印命令行参数
    // int i = 0;
    // printf("argc=%d\n",argc);
    // for(i=0;i<argc;i++)
    //     printf("argv[%d]:%s\n",i,argv[i]);
    H5W(argc,argv);
    
    return 0;
}