#include <petsc.h>
int main(int argc,char **argv)
{
    Vec            x1, x2, *x3ts, *x4ts;
    Vec            x1r, x2r, x3r, x4r;
    PetscViewer    viewer;
    PetscRandom    rand;
    PetscMPIInt    rank;
    PetscInt       i, n = 6, n_timesteps = 5;
    PetscBool      equal;
    MPI_Comm       comm;
    PetscErrorCode ierr;

    // Mat            A;
    // PetscInt       *dnnz,*onnz,rstart,rend,M,N;
    
    
    PetscFunctionBegin;
    ierr = PetscInitialize(&argc, &argv, (char*) 0, NULL);if (ierr) return ierr;
    comm = PETSC_COMM_WORLD;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL, "-n", &n, NULL);
    ierr = PetscOptionsGetInt(NULL,NULL, "-n_timesteps", &n_timesteps, NULL);
    if (n_timesteps < 0) SETERRQ(comm, PETSC_ERR_USER_INPUT, \
                            "-n_timesteps must be nonnegative");

    PetscRandomCreate(comm,&rand);
    PetscRandomSetFromOptions(rand);

    // 在一个h5文件中的顶级group中写入一个vec
    PetscViewerHDF5Open(comm,"test.h5",FILE_MODE_WRITE,&viewer);
    VecCreate(comm,&x1);
    PetscObjectSetName((PetscObject)x1,"x1");
    VecSetSizes(x1,PETSC_DECIDE,n);
    VecSetFromOptions(x1);
    VecSetRandom(x1,rand);
    VecView(x1,viewer);
    // 继续写入一个vec作为该目录下的另一个dataset
    PetscObjectSetName((PetscObject)x1,"x1.2");
    VecView(x1,viewer);
    // 继续写入一个vec作为该目录下的另一个dataset
    PetscObjectSetName((PetscObject)x1,"x1.3");
    VecView(x1,viewer);

    // 在一个h5文件中的一个group写入一个vec
    PetscViewerHDF5PushGroup(viewer,"/test1");
    VecCreate(comm,&x2);
    PetscObjectSetName((PetscObject)x2,"x2");
    VecSetSizes(x2,PETSC_DECIDE,n);
    VecSetBlockSize(x2,2);
    VecSetFromOptions(x2);
    VecSetRandom(x2,rand);
    VecView(x2,viewer);
    // 在该group下再写一个vec
    PetscObjectSetName((PetscObject)x2,"x2.1");
    VecView(x2,viewer);
    PetscViewerHDF5PopGroup(viewer);

    // 在一个h5文件中的另一个group写入一个vec
    PetscViewerHDF5PushGroup(viewer,"/test2");
    VecCreate(comm,&x2);
    PetscObjectSetName((PetscObject)x2,"x2");
    VecSetSizes(x2,PETSC_DECIDE,n);
    VecSetBlockSize(x2,2);
    VecSetFromOptions(x2);
    VecSetRandom(x2,rand);
    VecView(x2,viewer);
    PetscViewerHDF5PopGroup(viewer);
    
    PetscViewerDestroy(&viewer);
    PetscFinalize();
    return 0;

}
