#include <petsc.h>

int main(int argc,char **argv)
{
    Vec            g_b,g_t,g_l,g_r;
    Vec            h_b,h_t,h_l,h_r;
    Vec            gg;
    PetscInt       low,high,ldim,iglobal;
    PetscViewer    viewer;
    PetscRandom    random;
    PetscMPIInt    rank;
    PetscInt       i, n = 10, n_timesteps = 5;
    MPI_Comm       comm;
    PetscErrorCode ierr;
    PetscInt       block=1;
    PetscInt       choose,thres=50;
    PetscScalar    rvalue,zero=0;
    PetscChar      fname[PETSC_MAX_PATH_LEN]="test.hdf5";

    ierr = PetscInitialize(&argc, &argv, (char*) 0, NULL);if (ierr) return ierr;
    comm = PETSC_COMM_WORLD;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL, "-n", &n, NULL);
    ierr = PetscOptionsGetInt(NULL,NULL, "-n_timesteps", &n_timesteps, NULL);
    if (n_timesteps < 0) SETERRQ(comm, PETSC_ERR_USER_INPUT, "-n_timesteps must be nonnegative");
    ierr = PetscOptionsGetInt(NULL,NULL, "-block", &block, NULL);
    PetscOptionsGetString(NULL,NULL,"-fname",fname,sizeof(fname),NULL);

    PetscRandomCreate(comm,&random);
    PetscRandomSetFromOptions(random);

    VecCreate(comm,&g_b);
    VecSetSizes(g_b,PETSC_DECIDE,n+1);
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

    for (i=0;i<ldim;i++)
    {
        iglobal=i+low;
        choose = rand()%100;       //根据0~100的随机数是否大于50来确定赋值给u还是h
        rvalue = (rand()%98+1)/100.0; //赋的值为0.01~0.99
        if(choose>=thres)  
        {
            VecSetValues(g_b,1,&iglobal,&rvalue,INSERT_VALUES);
            VecSetValues(h_b,1,&iglobal,&zero,INSERT_VALUES);
        }
        else
        {
            VecSetValues(g_b,1,&iglobal,&zero,INSERT_VALUES);
            VecSetValues(h_b,1,&iglobal,&rvalue,INSERT_VALUES);
        }

    }

    for (i=0;i<ldim;i++)
    {
        iglobal=i+low;
        choose = rand()%100;       //根据0~100的随机数是否大于50来确定赋值给u还是h
        rvalue = (rand()%98+1)/100.0; //赋的值为0.01~0.99
        if(choose>=thres)  
        {
            VecSetValues(g_t,1,&iglobal,&rvalue,INSERT_VALUES);
            VecSetValues(h_t,1,&iglobal,&zero,INSERT_VALUES);
        }
        else
        {
            VecSetValues(g_t,1,&iglobal,&zero,INSERT_VALUES);
            VecSetValues(h_t,1,&iglobal,&rvalue,INSERT_VALUES);
        }

    }
    for (i=0;i<ldim;i++)
    {
        iglobal=i+low;
        choose = rand()%100;       //根据0~100的随机数是否大于50来确定赋值给u还是h
        rvalue = (rand()%98+1)/100.0; //赋的值为0.01~0.99
        if(choose>=thres)  
        {
            VecSetValues(g_l,1,&iglobal,&rvalue,INSERT_VALUES);
            VecSetValues(h_l,1,&iglobal,&zero,INSERT_VALUES);
        }
        else
        {
            VecSetValues(g_l,1,&iglobal,&zero,INSERT_VALUES);
            VecSetValues(h_l,1,&iglobal,&rvalue,INSERT_VALUES);
        }

    }
    for (i=0;i<ldim;i++)
    {
        iglobal=i+low;
        choose = rand()%100;       //根据0~100的随机数是否大于50来确定赋值给u还是h
        rvalue = (rand()%98+1)/100.0; //赋的值为0.01~0.99
        if(choose>=thres)  
        {
            VecSetValues(g_r,1,&iglobal,&rvalue,INSERT_VALUES);
            VecSetValues(h_r,1,&iglobal,&zero,INSERT_VALUES);
        }
        else
        {
            VecSetValues(g_r,1,&iglobal,&zero,INSERT_VALUES);
            VecSetValues(h_r,1,&iglobal,&rvalue,INSERT_VALUES);
        }

    }
    VecAssemblyBegin(g_b);
    VecAssemblyEnd(g_b);
    VecAssemblyBegin(g_t);
    VecAssemblyEnd(g_t);

    // 在一个h5文件中的一个group写入一个vec
    PetscViewerHDF5Open(PETSC_COMM_WORLD,fname,FILE_MODE_WRITE,&viewer);
    PetscViewerHDF5PushGroup(viewer,"/boundary");

    PetscObjectSetName((PetscObject)g_b,"g_bottom");
    //VecSetBlockSize(g_b,block);
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

    PetscViewerHDF5PushGroup(viewer,"/g_2D");
    VecCreate(comm,&gg);
    PetscObjectSetName((PetscObject)gg,"g_2DInit");
    VecSetSizes(gg,PETSC_DECIDE,(n+1)*(n+1));
   //VecSetBlockSize(gg,block);
    VecSetFromOptions(gg);
    VecSetRandom(gg,random);
    VecView(gg,viewer); 
    PetscViewerHDF5PopGroup(viewer);

    VecDestroy(&g_b);VecDestroy(&g_t);VecDestroy(&g_r);VecDestroy(&g_l);
    VecDestroy(&h_b);VecDestroy(&h_t);VecDestroy(&h_r);VecDestroy(&h_l);
    VecDestroy(&gg);
    PetscViewerDestroy(&viewer);

    PetscFinalize();
    return 0;

}
