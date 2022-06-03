static char help[] = "Generate the input data. \n \
            -bdry 1: all of boundary are given random u\n\
            -bdry 2: all of boundary are given random h\n\
            -bdry 3: all of boundary are given random u or h.\n" ;

#include <petsc.h>

int main(int argc,char **argv)
{
    Vec            g_b,g_t,g_l,g_r;
    Vec            h_b,h_t,h_l,h_r;
    Vec            u_0;
    PetscInt       low,high,ldim,iglobal;
    PetscViewer    viewer;
    PetscMPIInt    rank;
    PetscInt       i, n = 10;
    MPI_Comm       comm;
    PetscErrorCode ierr;
    PetscInt       bdry=1;
    PetscInt       choose,thres=50;
    PetscScalar    rand_g,rand_h,rand_u0,zero=0;
    PetscChar      fname[PETSC_MAX_PATH_LEN]="default.hdf5";

    ierr = PetscInitialize(&argc, &argv, (char*) 0, help);if (ierr) return ierr;
    comm = PETSC_COMM_WORLD;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL, "-n", &n, NULL);
    PetscOptionsGetInt(NULL,NULL,"-bdry",&bdry,NULL);
    PetscOptionsGetString(NULL,NULL,"-fname",fname,sizeof(fname),NULL);

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
       
    rand_u0 = (rand()%10000)/10.0;       //u0赋值为0~1000K

    VecGetOwnershipRange(g_b,&low,&high);
    VecGetLocalSize(g_b,&ldim); //获得向量在各进程的部分索引
    for (i=0;i<ldim;i++)
    {
        iglobal=i+low;
        rand_g = (rand()%10000)/10.0;        //g赋值为0~1000K
        rand_h = (rand()%1000)/100.0;        //h赋值为0~10
        if(bdry==1)
        {
            VecSetValues(g_b,1,&iglobal,&rand_g,INSERT_VALUES);
            VecSetValues(h_b,1,&iglobal,&zero,  INSERT_VALUES);
        }
        else if (bdry==2)
        {
            VecSetValues(h_b,1,&iglobal,&rand_h,INSERT_VALUES);
            VecSetValues(g_b,1,&iglobal,&zero,  INSERT_VALUES);            
        }
        else if (bdry==3)
        {
            choose=rand()%2;
            if(choose==0)
            {
                VecSetValues(g_b,1,&iglobal,&rand_g,INSERT_VALUES);
                VecSetValues(h_b,1,&iglobal,&zero  ,INSERT_VALUES);
            }
            else
            {
                VecSetValues(h_b,1,&iglobal,&rand_h,INSERT_VALUES);
                VecSetValues(g_b,1,&iglobal,&zero  ,INSERT_VALUES);                
            }
        }
    }

    for (i=0;i<ldim;i++)
    {
        iglobal=i+low;
        rand_g = (rand()%10000)/10.0;        //g赋值为0~1000K
        rand_h = (rand()%1000)/100.0;        //h赋值为0~10
        if(bdry==1)
        {
            VecSetValues(g_t,1,&iglobal,&rand_g,INSERT_VALUES);
            VecSetValues(h_t,1,&iglobal,&zero,  INSERT_VALUES);
        }
        else if (bdry==2)
        {
            VecSetValues(h_t,1,&iglobal,&rand_h,INSERT_VALUES);
            VecSetValues(g_t,1,&iglobal,&zero,  INSERT_VALUES);            
        }
        else if (bdry==3)
        {
            choose=rand()%2;
            if(choose==0)
            {
                VecSetValues(g_t,1,&iglobal,&rand_g,INSERT_VALUES);
                VecSetValues(h_t,1,&iglobal,&zero  ,INSERT_VALUES);
            }
            else
            {
                VecSetValues(h_t,1,&iglobal,&rand_h,INSERT_VALUES);
                VecSetValues(g_t,1,&iglobal,&zero  ,INSERT_VALUES);                
            }
        }
    }

    for (i=0;i<ldim;i++)
    {
        iglobal=i+low;
        rand_g = (rand()%10000)/10.0;        //g赋值为0~1000K
        rand_h = (rand()%1000)/100.0;        //h赋值为0~10
        if(bdry==1)
        {
            VecSetValues(g_r,1,&iglobal,&rand_g,INSERT_VALUES);
            VecSetValues(h_r,1,&iglobal,&zero,  INSERT_VALUES);
        }
        else if (bdry==2)
        {
            VecSetValues(h_r,1,&iglobal,&rand_h,INSERT_VALUES);
            VecSetValues(g_r,1,&iglobal,&zero,  INSERT_VALUES);            
        }
        else if (bdry==3)
        {
            choose=rand()%2;
            if(choose==0)
            {
                VecSetValues(g_r,1,&iglobal,&rand_g,INSERT_VALUES);
                VecSetValues(h_r,1,&iglobal,&zero  ,INSERT_VALUES);
            }
            else
            {
                VecSetValues(h_r,1,&iglobal,&rand_h,INSERT_VALUES);
                VecSetValues(g_r,1,&iglobal,&zero  ,INSERT_VALUES);                
            }
        }
    }

    for (i=0;i<ldim;i++)
    {
        iglobal=i+low;
        rand_g = (rand()%10000)/10.0;        //g赋值为0~1000K
        rand_h = (rand()%1000)/100.0;        //h赋值为0~10
        if(bdry==1)
        {
            VecSetValues(g_l,1,&iglobal,&rand_g,INSERT_VALUES);
            VecSetValues(h_l,1,&iglobal,&zero,  INSERT_VALUES);
        }
        else if (bdry==2)
        {
            VecSetValues(h_l,1,&iglobal,&rand_h,INSERT_VALUES);
            VecSetValues(g_l,1,&iglobal,&zero,  INSERT_VALUES);            
        }
        else if (bdry==3)
        {
            choose=rand()%2;
            if(choose==0)
            {
                VecSetValues(g_l,1,&iglobal,&rand_g,INSERT_VALUES);
                VecSetValues(h_l,1,&iglobal,&zero  ,INSERT_VALUES);
            }
            else
            {
                VecSetValues(h_l,1,&iglobal,&rand_h,INSERT_VALUES);
                VecSetValues(g_l,1,&iglobal,&zero  ,INSERT_VALUES);                
            }
        }
    }

    // VecAssemblyBegin(g_b);
    // VecAssemblyEnd(g_b);
    // VecAssemblyBegin(g_t);
    // VecAssemblyEnd(g_t);


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

    PetscViewerHDF5PushGroup(viewer,"/u0_2D");
    VecCreate(comm,&u_0);
    PetscObjectSetName((PetscObject)u_0,"u0_2DInit");
    VecSetSizes(u_0,PETSC_DECIDE,n*n);
    VecSetFromOptions(u_0);
    VecGetOwnershipRange(u_0,&low,&high);
    VecGetLocalSize(u_0,&ldim); 
    for (i=0;i<ldim;i++)
    {
        iglobal=i+low;
        rand_u0 = (rand()%10000)/10.0;   
        VecSetValues(u_0,1,&iglobal,&rand_u0,INSERT_VALUES);
    }
   //VecSetBlockSize(gg,block);
    VecView(u_0,viewer); 
    PetscViewerHDF5PopGroup(viewer);

    VecDestroy(&g_b);VecDestroy(&g_t);VecDestroy(&g_r);VecDestroy(&g_l);
    VecDestroy(&h_b);VecDestroy(&h_t);VecDestroy(&h_r);VecDestroy(&h_l);
    VecDestroy(&u_0);
    PetscViewerDestroy(&viewer);

    PetscFinalize();
    return 0;

}
