static char help[] = "Generate the input data. \n \
            -bdry 1: all of boundary are given u\n\
            -bdry 2: all of boundary are given h\n\
            -bdry 3: all of boundary are given u or h.\n";

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
    PetscInt       choose,thres=50;
    PetscScalar    rand_u,rand_h,zero=0;
    PetscChar      fname[PETSC_MAX_PATH_LEN]="default.hdf5";

    ierr = PetscInitialize(&argc, &argv, (char*) 0, help);if (ierr) return ierr;
    comm = PETSC_COMM_WORLD;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL, "-n", &n, NULL);
    PetscOptionsGetInt(NULL,NULL,"-bdry",&bdry,NULL);
    PetscOptionsGetString(NULL,NULL,"-fname",fname,sizeof(fname),NULL);

    PetscRandomCreate(comm,&random);
    PetscRandomSetFromOptions(random);

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
    for (i=0;i<ldim;i++)
    {
        rand_u = rand()%273;        //u赋值为0~273
        rand_h = rand()%1000/100.0; //h赋值为0~10
        iglobal=i+low;
        if(bdry==1)
        {
            VecSetValues(g_b,1,&iglobal,&rand_u,INSERT_VALUES);
        }
        choose = rand()%100;        //根据0~100的随机数是否大于50来确定赋值给u还是h
        rand_u = rand()%273;        //u赋值为0~273
        rand_h = rand()%1000/100.0; //h赋值为0~10
        if(choose>=thres)  
        {
            VecSetValues(g_b,1,&iglobal,&rand_u,INSERT_VALUES);
            VecSetValues(h_b,1,&iglobal,&zero,INSERT_VALUES);
        }
        else
        {
            VecSetValues(g_b,1,&iglobal,&zero,INSERT_VALUES);
            VecSetValues(h_b,1,&iglobal,&rand_h,INSERT_VALUES);
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
