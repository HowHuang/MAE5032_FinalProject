static char help[] = "Generate the fixed input data. \n \
            -bdry 1: all of boundary are given invariant g\n\
            -bdry 2: all of boundary are given invariant h.\n \
            -gv: specify the invariant value g\n \
            -hv: specify the invariant value h\n \
            -u0: specify the init(T=0) value of u ";

#include <petsc.h>

int main(int argc,char **argv)
{
    Vec            g_b,g_t,g_l,g_r;
    Vec            h_b,h_t,h_l,h_r;
    Vec            u_0;
    PetscInt       low,high,ldim,iglobal;
    PetscViewer    viewer;
    PetscRandom    random;
    PetscMPIInt    rank;
    PetscInt       i, n = 10;
    MPI_Comm       comm;
    PetscErrorCode ierr;
    PetscInt       bdry=1;
    PetscScalar    g_value,h_value,u0_value,zero=0;
    PetscChar      fname[PETSC_MAX_PATH_LEN]="default.hdf5";


    ierr = PetscInitialize(&argc, &argv, (char*) 0, help);if (ierr) return ierr;
    comm = PETSC_COMM_WORLD;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);

    u0_value = (rand()%10000)/10.0;       //u0赋值为0~1000K
    g_value = (rand()%10000)/10.0;        //g赋值为0~1000K
    h_value = (rand()%1000)/100.0;        //h赋值为0~10
    PetscPrintf(comm,"u0value: %g, g_value: %g, h_value: %g\n",u0_value,g_value,h_value);

    PetscOptionsGetInt(NULL,NULL, "-n", &n, NULL);
    PetscOptionsGetInt(NULL,NULL,"-bdry",&bdry,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-gv",&g_value,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-hv",&h_value,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-u0",&u0_value,NULL);
    PetscOptionsGetString(NULL,NULL,"-fname",fname,sizeof(fname),NULL);

    if(bdry==1)
    {
        PetscPrintf(comm,"use the boundary g : %g\n",g_value);
    }
    else if(bdry==2)
    {
        PetscPrintf(comm,"use the boundary h : %g\n",h_value);
    }
    else
    {
        PetscPrintf(comm,"Please specify boundray type: 1 for u or 2 for h\n");
        return -1;
    }
    PetscPrintf(comm,"use the init u0 : %g\n",u0_value);

    VecCreate(comm,&g_b);
    VecSetSizes(g_b,PETSC_DECIDE,n);
    VecSetFromOptions(g_b);
    VecDuplicate(g_b,&h_b);

    VecDuplicate(g_b,&g_t);
    VecDuplicate(g_b,&h_t);
    VecDuplicate(g_b,&g_l);
    VecDuplicate(g_b,&h_l);
    VecDuplicate(g_b,&g_r);
    VecDuplicate(g_b,&h_r);

    VecGetOwnershipRange(g_b,&low,&high);
    VecGetLocalSize(g_b,&ldim); //获得向量在各进程的部分索引

    if(bdry==1)
    {
        VecSet(g_b,g_value);
        VecSet(g_t,g_value);
        VecSet(g_l,g_value);
        VecSet(g_r,g_value);

        VecSet(h_b,zero);
        VecSet(h_t,zero);
        VecSet(h_l,zero);
        VecSet(h_r,zero);
    }
    else if (bdry==2)
    {
        VecSet(h_b,h_value);
        VecSet(h_t,h_value);
        VecSet(h_l,h_value);
        VecSet(h_r,h_value);

        VecSet(g_b,zero);
        VecSet(g_t,zero);
        VecSet(g_l,zero);
        VecSet(g_r,zero);        
    }

    VecAssemblyBegin(g_b);
    VecAssemblyEnd(g_b);
    VecAssemblyBegin(g_t);
    VecAssemblyEnd(g_t);

    PetscViewerHDF5Open(PETSC_COMM_WORLD,fname,FILE_MODE_WRITE,&viewer);
    PetscViewerHDF5PushGroup(viewer,"/boundary");

    PetscObjectSetName((PetscObject)g_b,"g_bottom");
    VecView(g_b,viewer);
    PetscObjectSetName((PetscObject)h_b,"h_bottom");
    VecView(h_b,viewer);

    PetscObjectSetName((PetscObject)g_t,"g_top");
    VecView(g_t,viewer);   
    PetscObjectSetName((PetscObject)h_t,"h_top");
    VecView(h_t,viewer);  

    PetscObjectSetName((PetscObject)g_l,"g_left");
    VecView(g_l,viewer);  
    PetscObjectSetName((PetscObject)h_l,"h_left");
    VecView(h_l,viewer);  

    PetscObjectSetName((PetscObject)g_r,"g_right");
    VecView(g_r,viewer);  
    PetscObjectSetName((PetscObject)h_r,"h_right");
    VecView(h_r,viewer);  
    PetscViewerHDF5PopGroup(viewer);

    PetscViewerHDF5PushGroup(viewer,"/u0_2D");
    VecCreate(comm,&u_0);
    PetscObjectSetName((PetscObject)u_0,"u0_2DInit");
    VecSetSizes(u_0,PETSC_DECIDE,n*n);
    VecSetFromOptions(u_0);
    VecSet(u_0,u0_value);
    VecView(u_0,viewer); 
    PetscViewerHDF5PopGroup(viewer);

    VecDestroy(&g_b);VecDestroy(&g_t);VecDestroy(&g_r);VecDestroy(&g_l);
    VecDestroy(&h_b);VecDestroy(&h_t);VecDestroy(&h_r);VecDestroy(&h_l);
    VecDestroy(&u_0);
    PetscViewerDestroy(&viewer);

    PetscFinalize();
    return 0;

}
