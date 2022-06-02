#include <petsc.h>

int main(int argc,char **argv)
{
    Vec            g_b,g_t,g_l,g_r;
    Vec            h_b,h_t,h_l,h_r;
    Vec            gg;
    PetscInt       low,high,ldim,iglobal,lsize;
    PetscViewer    viewer;
    PetscRandom    random;
    PetscMPIInt    rank,size;
    PetscInt       i, n = 10;
    MPI_Comm       comm;
    PetscErrorCode ierr;
    PetscInt       block=1;
    PetscInt       choose,thres=50;
    PetscScalar    rvalue,zero=0;
    PetscChar      fname[PETSC_MAX_PATH_LEN]="test.hdf5";
    PetscBool      vstage2=PETSC_FALSE;
    PetscBool      vstage3=PETSC_FALSE;

    ierr = PetscInitialize(&argc, &argv, (char*) 0, NULL);if (ierr) return ierr;
    comm = PETSC_COMM_WORLD;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
    ierr = MPI_Comm_size(comm, &size);CHKERRMPI(ierr);

    ierr = PetscOptionsGetInt(NULL,NULL, "-n", &n, NULL);
    ierr = PetscOptionsGetInt(NULL,NULL, "-block", &block, NULL);
    PetscOptionsGetString(NULL,NULL,"-fname",fname,sizeof(fname),NULL);
    PetscOptionsGetBool(NULL,NULL,"-sizes_set",&vstage2,NULL);
    PetscOptionsGetBool(NULL,NULL,"-type_set",&vstage3,NULL);

    PetscViewerHDF5Open(PETSC_COMM_WORLD,fname,FILE_MODE_READ,&viewer);
    PetscViewerHDF5PushGroup(viewer,"/boundary");

    VecCreate(comm,&g_b);
    PetscObjectSetName((PetscObject)g_b, "g_bottom");
    VecCreate(comm,&h_b);
    PetscObjectSetName((PetscObject)h_b, "h_bottom");
    // if (vstage2) 
    // {
    //     PetscPrintf(PETSC_COMM_WORLD,"Setting vector sizes...\n");
    //     if (size > 1) 
    //     {
    //         if (rank == 0) 
    //         {
    //             lsize = (n+1)/size + size;
    //             VecSetSizes(g_b,lsize,n+1);
    //         } 
    //         else if (rank == size-1) 
    //         {
    //             lsize = (n+1)/size - size;
    //             VecSetSizes(g_b,lsize,n+1);
    //         } 
    //         else 
    //         {
    //             lsize = (n+1)/size;
    //             VecSetSizes(g_b,lsize,n+1);
    //         }
    //     } 
    //     else 
    //     {
    //         VecSetSizes(g_b,n+1,n+1);
    //     }
    // }

    // if (vstage3) 
    // {
    //     PetscPrintf(PETSC_COMM_WORLD,"Setting vector type...\n");
    //     VecSetType(g_b, VECMPI);
    // }   

    VecLoad(g_b,viewer);
    VecLoad(h_b,viewer);
    PetscViewerHDF5PopGroup(viewer);
    PetscViewerDestroy(&viewer);
    VecView(g_b,PETSC_VIEWER_STDOUT_WORLD);

    // VecDestroy(&g_b);VecDestroy(&g_t);VecDestroy(&g_r);VecDestroy(&g_l);
    // VecDestroy(&h_b);VecDestroy(&h_t);VecDestroy(&h_r);VecDestroy(&h_l);
    // VecDestroy(&gg);
    // PetscViewerDestroy(&viewer);

    PetscFinalize();
    return 0;

}
