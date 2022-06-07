#include <petsc.h>

void H5W(int argc, char** argv)
{
    PetscMPIInt         rank, size;
    PetscInt            i,m=20,low,high,ldim,iglobal,lsize;
    PetscScalar         v;
    Vec                 u;
    PetscViewer         viewer;
    PetscBool           ishdf5 = PETSC_TRUE;
    PetscScalar const*  values;
    PetscChar           fname[PETSC_MAX_PATH_LEN]="default.hdf5";
    PetscChar           oname[PETSC_MAX_PATH_LEN]="default";
    #if defined(PETSC_USE_LOG)
    PetscLogEvent       VECTOR_GENERATE,VECTOR_READ;
    #endif
    
    PetscInitialize(&argc, &argv, (char*)0, NULL);
    PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);
    PetscOptionsGetString(NULL,NULL,"-fname",fname,sizeof(fname),NULL);
    PetscOptionsGetString(NULL,NULL,"-oname",oname,sizeof(oname),NULL);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    // Generate vector, then write it in the given data format 
    PetscLogEventRegister("Generate Vector",VEC_CLASSID,&VECTOR_GENERATE);
    PetscLogEventBegin(VECTOR_GENERATE,0,0,0,0);

    VecCreate(PETSC_COMM_WORLD,&u); //创建
    PetscObjectSetName((PetscObject)u, oname); //起名字
    VecSetSizes(u,PETSC_DECIDE,m);   //设置全局大小，局部由petsc决定
    VecSetFromOptions(u);   //设置为可以由命令行指定
    VecGetOwnershipRange(u,&low,&high); //获得向量在各进程的部分索引
    VecGetLocalSize(u,&ldim);   //获得本地size

    for(i=0;i<ldim;i++)     //遍历本地
    {
        iglobal = i+low;    //获得全局索引
        v = (PetscScalar)(i+low);   //值等于索引
        VecSetValues(u,1,&iglobal,&v,INSERT_VALUES);    //赋值
    }   

    VecAssemblyBegin(u);
    VecAssemblyEnd(u);
    VecView(u,PETSC_VIEWER_STDOUT_WORLD);   //查看

    // 把vec输出到hdf5文件中，FILE_MODE_WRITE表示如果文件不存在则创建
    PetscPrintf(PETSC_COMM_WORLD,"writing vector in hdf5 to vector.dat ...\n");
    PetscViewerHDF5Open(PETSC_COMM_WORLD,fname,FILE_MODE_WRITE,&viewer);
    VecView(u,viewer);
    PetscViewerDestroy(&viewer);
    VecDestroy(&u);

    PetscLogEventEnd(VECTOR_GENERATE,0,0,0,0);
}