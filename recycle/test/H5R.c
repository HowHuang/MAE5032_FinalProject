#include <petsc.h>

void H5R(int argc, char** argv)
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
    PetscBool           vstage2=PETSC_FALSE;
    PetscBool           vstage3=PETSC_FALSE;
    
    PetscInitialize(&argc, &argv, (char*)0, NULL);
    PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);
    PetscOptionsGetString(NULL,NULL,"-fname",fname,sizeof(fname),NULL);
    PetscOptionsGetString(NULL,NULL,"-oname",oname,sizeof(oname),NULL);
    PetscOptionsGetBool(NULL,NULL,"-sizes_set",&vstage2,NULL);
    PetscOptionsGetBool(NULL,NULL,"-type_set",&vstage3,NULL);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    // Generate vector, then write it in the given data format 
    PetscLogEventRegister("Read Vector",VEC_CLASSID,&VECTOR_READ);
    PetscLogEventBegin(VECTOR_READ,0,0,0,0);
    PetscPrintf(PETSC_COMM_WORLD,"reading vector in hdf5 from vector.dat ...\n");
    PetscViewerHDF5Open(PETSC_COMM_WORLD,fname,FILE_MODE_READ,&viewer);

    VecCreate(PETSC_COMM_WORLD,&u);
    PetscObjectSetName((PetscObject)u, oname);

    if (vstage2) 
    {
        PetscPrintf(PETSC_COMM_WORLD,"Setting vector sizes...\n");
        if (size > 1) 
        {
            if (rank == 0) 
            {
                lsize = m/size + size;
                VecSetSizes(u,lsize,m);
            } 
            else if (rank == size-1) 
            {
                lsize = m/size - size;
                VecSetSizes(u,lsize,m);
            } 
            else 
            {
                lsize = m/size;
                VecSetSizes(u,lsize,m);
            }
        } 
        else 
        {
            VecSetSizes(u,m,m);
        }
    }

    if (vstage3) 
    {
        PetscPrintf(PETSC_COMM_WORLD,"Setting vector type...\n");
        VecSetType(u, VECMPI);
    }   

    VecLoad(u,viewer);
    PetscViewerDestroy(&viewer);
    PetscLogEventEnd(VECTOR_READ,0,0,0,0);
    VecView(u,PETSC_VIEWER_STDOUT_WORLD);

    VecGetArrayRead(u,&values);
    VecGetLocalSize(u,&ldim);
    VecGetOwnershipRange(u,&low,NULL);
    for (i=0;i<ldim;i++)
    {

    }
    VecRestoreArrayRead(u,&values);

    VecDestroy(&u);
    PetscFinalize();
    return 0;

}