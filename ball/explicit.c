#include <petsc.h>
#include "petscviewerhdf5.h" 
#include <hdf5.h>
#include "plicit.h"

// 根据给定的边界条件和物理参数，计算对应点关联系数矩阵行及右手向量中的元素值
void ExplicitIterationMaterial(InputPara* IP, Bound* bound, IterMaterial* IM, enum Location loc)
{
    //~ InputPara       Known;
    //~ bound           Known;
    //~ IterMaterial    Unknown;
    // 计算右手向量b中的温度热源部分、热流热源部分和热补给部分。
    PetscScalar base = (IP->k * IP->dt) / (IP->rho * IP->c * PetscPowScalar(IP->dl,2));
    PetscScalar partF = IP->f * IP->dt / (IP->rho * IP->c);
    PetscScalar partH = bound->hb * IP->dt / (IP->rho * IP->c * pow(IP->dl,2)) + \
                        bound->ht * IP->dt / (IP->rho * IP->c * pow(IP->dl,2)) + \
                        bound->hl * IP->dt / (IP->rho * IP->c * pow(IP->dl,2)) + \
                        bound->hr * IP->dt / (IP->rho * IP->c * pow(IP->dl,2));
    PetscScalar partU = 2 * base * (bound->ut + bound->ur + bound->ub + bound->ul);      

    // 根据不同的位置对系数矩阵和右手向量进行特别计算              
    switch (loc)
    {
    case RightTop:
    {
        IM->W = base;
        IM->S = base;
        IM->b = partF + partU + partH;     

        if(bound->tr==1 && bound->tt==1)
            IM->P = 1-(1+1+2+2)*base;
        else if (bound->tr==2 && bound->tt==2)
            IM->P = 1-(1+1)*base;
        else 
            IM->P = 1-(1+1+2)*base;
    }
    break;

    case LeftTop:
    {
        IM->E = base;
        IM->S = base;
        IM->b = partF + partU + partH;   

        if(bound->tl==1 && bound->tt==1)
            IM->P = 1-(1+1+2+2)*base;
        else if (bound->tl==2 && bound->tt==2)
            IM->P = 1-(1+1)*base;
        else 
            IM->P = 1-(1+1+2)*base;         
    }
    break;

    case LeftBottom:
    {
        IM->E = base;
        IM->S = base;
        IM->b = partF + partU + partH;   
        if(bound->tl==1 && bound->tb==1)
            IM->P = 1-(1+1+2+2)*base;
        else if (bound->tl==2 && bound->tb==2)
            IM->P = 1-(1+1)*base;
        else 
            IM->P = 1-(1+1+2)*base;
    }
    break;

    case RightBottom:
    {
        IM->E = base;
        IM->S = base;
        IM->b = partF + partU + partH;  

        if(bound->tr==1 && bound->tb==1)
            IM->P = 1-(1+1+2+2)*base;
        else if (bound->tr==2 && bound->tb==2)
            IM->P = 1-(1+1)*base;
        else 
            IM->P = 1-(1+1+2)*base;    
    }
    break;

    case RightSide:
    {
        IM->W = base;
        IM->S = base;
        IM->N = base;
        IM->b = partF + partU + partH;   

        if(bound->tr==1)
            IM->P = 1-(1+1+1+2)*base;    
        else
            IM->P = 1-(1+1+1)*base;   
    }
    break;

    case TopSide:
    {
        IM->W = base;
        IM->E = base;
        IM->S = base;
        IM->b = partF + partU + partH;   

        if(bound->tt==1)
            IM->P = 1-(1+1+1+2)*base;    
        else
            IM->P = 1-(1+1+1)*base;   
    }
    break;

    case LeftSide:
    {
        IM->N = base;
        IM->E = base;
        IM->S = base;
        IM->b = partF + partU + partH;  

        if(bound->tl==1)
            IM->P = 1-(1+1+1+2)*base;    
        else
            IM->P = 1-(1+1+1)*base;    
    }   
    break;

    case BottomSide:
    {
        IM->N = base;
        IM->E = base;
        IM->W = base;
        IM->b = partF + partU + partH;   

        if(bound->tb==1)
            IM->P = 1-(1+1+1+2)*base;    
        else
            IM->P = 1-(1+1+1)*base; 
    }    
    break;

    case Internal:
    {
        IM->b = partF;   
        IM->N = base;
        IM->E = base;
        IM->W = base;
        IM->S = base;
        IM->P = 1-4*base;
    }    
    break;

    default:
        printf("Please specify the location\n");
    }

}
    
int Explicit(int argc,char **argv)
{
    Vec             u_0,b,u_t,u_tplus,temp;     //DIM = (n*n) x 1
    Mat             A;                          //DIM = (n*n) x (n*n)
    PetscViewer     viewer;    
    PetscMPIInt     rank,size;
    PetscInt        i, j, r, n = 10;
    PetscInt        restart;
    MPI_Comm        comm;
    PetscErrorCode  ierr;
    PetscChar       fname[PETSC_MAX_PATH_LEN]="input.hdf5";
    PetscChar       dsname[PETSC_MAX_PATH_LEN]="default";
    PetscInt        col[5];
    PetscScalar     value[5];

    InputPara       IP; 
    Bound           bound;
    IterMaterial    IM;
    Location        loc;

    PetscScalar     *g_b, *g_t, *g_l, *g_r;
    PetscScalar     *h_b, *h_t, *h_l, *h_r;
    PetscScalar     *t_b, *t_t, *t_l, *t_r;

    PetscScalar     paras[7];
    
    // iteration parts
    PetscInt        its=0 , maxIts = 100, maxItsW = 1000;
    Vec             *u; 
    Vec             step;
    PetscInt        istep;        

    ierr = PetscInitialize(&argc, &argv, (char*) 0, NULL);if (ierr) return ierr;
    comm = PETSC_COMM_WORLD;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
    ierr = MPI_Comm_size(comm, &size);CHKERRMPI(ierr);

    ierr = PetscOptionsGetString(NULL,NULL,"-fname",fname,sizeof(fname),NULL);
    ierr = PetscOptionsGetInt(NULL,NULL, "-maxIts", &maxIts, NULL);
    ierr = PetscOptionsGetInt(NULL,NULL, "-maxItsW", &maxItsW, NULL);
    ierr = PetscOptionsGetInt(NULL,NULL, "-restart", &restart, NULL);

    hid_t   file_id, group_id, dataset_id;

    file_id=H5Fopen(fname,H5F_ACC_RDONLY,H5P_DEFAULT);

    group_id=H5Gopen(file_id,"/Parameters",H5P_DEFAULT);
    dataset_id=H5Dopen(group_id,"paras",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,paras);
    H5Dclose(dataset_id);
    H5Gclose(group_id);

    IP.dt=paras[0]; IP.dl=paras[1]; IP.rho=paras[2];
    IP.c=paras[3];  IP.k=paras[4];  IP.f=paras[5];

    n = (PetscInt)(paras[6]);

    PetscMalloc4(n, &g_b, n, &g_t, n, &g_l, n, &g_r);
    PetscMalloc4(n, &h_b, n, &h_t, n, &h_l, n, &h_r);
    PetscMalloc4(n, &t_b, n, &t_t, n, &t_l, n, &t_r);

    group_id=H5Gopen(file_id,"/boundary",H5P_DEFAULT);

    dataset_id=H5Dopen(group_id,"g_bottom",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,g_b);
    H5Dclose(dataset_id);
    dataset_id=H5Dopen(group_id,"g_top",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,g_t);
    H5Dclose(dataset_id);
    dataset_id=H5Dopen(group_id,"g_left",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,g_l);
    H5Dclose(dataset_id);
    dataset_id=H5Dopen(group_id,"g_right",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,g_r);
    H5Dclose(dataset_id);

    dataset_id=H5Dopen(group_id,"h_bottom",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,h_b);
    H5Dclose(dataset_id);
    dataset_id=H5Dopen(group_id,"h_top",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,h_t);
    H5Dclose(dataset_id);
    dataset_id=H5Dopen(group_id,"h_left",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,h_l);
    H5Dclose(dataset_id);
    dataset_id=H5Dopen(group_id,"h_right",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,h_r);
    H5Dclose(dataset_id);

    dataset_id=H5Dopen(group_id,"t_bottom",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,t_b);
    H5Dclose(dataset_id);
    dataset_id=H5Dopen(group_id,"t_top",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,t_t);
    H5Dclose(dataset_id);
    dataset_id=H5Dopen(group_id,"t_left",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,t_l);
    H5Dclose(dataset_id);
    dataset_id=H5Dopen(group_id,"t_right",H5P_DEFAULT);
    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,t_r);
    H5Dclose(dataset_id);

    H5Gclose(group_id);

    H5Fclose(file_id);

    // ~ Load data from specified hdf5 file.
    PetscViewerHDF5Open(PETSC_COMM_WORLD,fname,FILE_MODE_READ,&viewer);
    
    PetscViewerHDF5PushGroup(viewer,"/Init");
    VecCreate(comm,&u_0);
    PetscObjectSetName((PetscObject)u_0, "u_0");
    VecLoad(u_0,viewer);
    PetscViewerHDF5PopGroup(viewer);

    // ~ Load finished

    // ~ Generate sparse matrix A (coefficient matrix)
    // ~ A is a pentadiagonal matrix
    MatCreate(comm,&A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n*n,n*n);
    MatSetFromOptions(A);
    MatSetUp(A);
    
    // ~ Generate the vector b for Au_t+b=u_tplus
    VecCreate(comm,&b);
    VecSetSizes(b,PETSC_DECIDE,n*n);
    VecSetFromOptions(b);
    VecSetUp(b);

    VecDuplicate(b,&u_t);
    VecDuplicate(b,&u_tplus);
    VecDuplicate(b,&temp);
    VecDuplicateVecs(b,maxIts+1,&u);  // store u for each timestep

    // ~ Assign values to A and b
 
    if(rank==0)
    {
        for(i=0;i<n;i++)            //global row index of u (consider as matrix)        
        {   
            for(j=0;j<n;j++)        //global col index of u (consider as matrix) 
            {
                r=i*n+j;            //global row index of A (same as b)
                if(i==0)
                {
                    if(j==0)            //! 左上角的点
                    {   
                        loc = LeftTop;
                        bound.ut=g_t[j];bound.ht=h_t[j];bound.tt=t_t[j];
                        bound.ul=g_l[i];bound.hl=h_l[i];bound.tl=t_l[i];
                        bound.ub=0;     bound.hb=0;     bound.tb=0;
                        bound.ur=0;     bound.hr=0;     bound.tr=0;

                        ExplicitIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r;col[1]=r+1;col[2]=r+n;
                        value[0]=IM.P;value[1]=IM.E;value[2]=IM.S;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES); 
                        VecSetValue(b,r,IM.b,INSERT_VALUES);
                    }
                    else if(j==n-1)     //! 右上角的点
                    {
                        loc = RightTop;
                        bound.ut=g_t[j];bound.ht=h_t[j];bound.tt=t_t[j];
                        bound.ur=g_r[i];bound.hr=h_r[i];bound.tr=t_r[i];
                        bound.ub=0;     bound.hb=0;     bound.tb=0;
                        bound.ul=0;     bound.hl=0;     bound.tl=0;
                        ExplicitIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-1;col[1]=r;col[2]=r+n;
                        value[0]=IM.W;value[1]=IM.P;value[2]=IM.S;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);   
                        VecSetValue(b,r,IM.b,INSERT_VALUES);       
                    }
                    else                //! 上边界的点（除顶点外）
                    {   
                        loc = TopSide;
                        bound.ut=g_t[j];bound.ht=h_t[j];bound.tt=t_t[j];
                        bound.ur=0;     bound.hr=0;     bound.tr=0;
                        bound.ub=0;     bound.hb=0;     bound.tb=0;
                        bound.ul=0;     bound.hl=0;     bound.tl=0;
                        ExplicitIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-1;col[1]=r;col[2]=r+1;col[3]=r+n;
                        value[0]=IM.W;value[1]=IM.P;value[2]=IM.E;value[3]=IM.S;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);
                        VecSetValue(b,r,IM.b,INSERT_VALUES);              
                    }
                }
                else if(i==n-1)
                {
                    if(j==0)            //! 左下角的点
                    {
                        loc = LeftBottom;
                        bound.ub=g_b[j];bound.hb=h_b[j];bound.tb=t_b[j];
                        bound.ur=0;     bound.hr=0;     bound.tr=0;
                        bound.ut=0;     bound.ht=0;     bound.tt=0;
                        bound.ul=g_l[i];bound.hl=h_l[i];bound.tl=t_l[i];
                        ExplicitIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-n;col[1]=r;col[2]=r+1;
                        value[0]=IM.N;value[1]=IM.P;value[2]=IM.E;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);  
                        VecSetValue(b,r,IM.b,INSERT_VALUES);                                   
                    }
                    else if(j==n-1)     //! 右下角的点
                    {   
                        loc = RightBottom;
                        bound.ub=g_b[j];bound.hb=h_b[j];bound.tb=t_b[j];
                        bound.ul=0;     bound.hl=0;     bound.tl=0;
                        bound.ut=0;     bound.ht=0;     bound.tt=0;
                        bound.ur=g_r[i];bound.hr=h_r[i];bound.tr=t_r[i];
                        ExplicitIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-n;col[1]=r-1;col[2]=r;
                        value[0]=IM.N;value[1]=IM.W;value[2]=IM.P;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);      
                        VecSetValue(b,r,IM.b,INSERT_VALUES);                                     
                    }
                    else                //! 下边界的点 （除顶点外）
                    {
                        loc = BottomSide;
                        bound.ub=g_b[j];bound.hb=h_b[j];bound.tb=t_b[j];
                        bound.ul=0;     bound.hl=0;     bound.tl=0;
                        bound.ut=0;     bound.ht=0;     bound.tt=0;
                        bound.ur=0;     bound.hr=0;     bound.tr=0;
                        ExplicitIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-n;col[1]=r-1;col[2]=r;col[3]=r+1;
                        value[0]=IM.N;value[1]=IM.W;value[2]=IM.P;value[3]=IM.E;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);  
                        VecSetValue(b,r,IM.b,INSERT_VALUES);                                       
                    }
                }
                else
                {
                    if(j==0)            //! 左边界的点（除顶点外）
                    {
                        loc = LeftSide;
                        bound.ub=0;     bound.hb=0;     bound.tb=0;
                        bound.ul=g_l[i];bound.hl=h_l[i];bound.tl=t_l[i];
                        bound.ut=0;     bound.ht=0;     bound.tt=0;
                        bound.ur=0;     bound.hr=0;     bound.tr=0;
                        ExplicitIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-n;col[1]=r;col[2]=r+1;col[3]=r+n;
                        value[0]=IM.N;value[1]=IM.P;value[2]=IM.E;value[3]=IM.S;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);
                        VecSetValue(b,r,IM.b,INSERT_VALUES);

                    }
                    else if(j==n-1)     //! 右边界的点（除顶点外）
                    {
                        loc = RightSide;
                        bound.ub=0;     bound.hb=0;     bound.tb=0;
                        bound.ur=g_r[i];bound.hr=h_r[i];bound.tr=t_r[i];
                        bound.ut=0;     bound.ht=0;     bound.tt=0;
                        bound.ul=0;     bound.hl=0;     bound.tl=0;
                        ExplicitIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-n;col[1]=r-1;col[2]=r;col[3]=r+n;
                        value[0]=IM.N;value[1]=IM.W;value[2]=IM.P;value[3]=IM.S;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);
                        VecSetValue(b,r,IM.b,INSERT_VALUES);                  
                    }
                    else                //! 内部点
                    {
                        loc = Internal;
                        bound.ub=0;     bound.hb=0;     bound.tb=0;
                        bound.ur=0;     bound.hr=0;     bound.tr=0;
                        bound.ut=0;     bound.ht=0;     bound.tt=0;
                        bound.ul=0;     bound.hl=0;     bound.tl=0;
                        ExplicitIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-n;col[1]=r-1;col[2]=r;col[3]=r+1;col[4]=r+n;
                        value[0]=IM.N;value[1]=IM.W;value[2]=IM.P;value[3]=IM.E;value[4]=IM.S;
                        MatSetValues(A,1,&r,5,col,value,INSERT_VALUES);
                        VecSetValue(b,r,IM.b,INSERT_VALUES);
                    }
                }
            }
        }
    }

    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); 
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    // ~ ---------------------------------------------------------------------------
    // ~ Now, we have A, u_t(u_0), b and u_tplus(unknown) Let's begin the iteration.
    // ~ We will store u_t at any timestep in a HDF5 file.
    VecCreate(comm,&step);
    VecSetSizes(step,PETSC_DECIDE,1);
    VecSetFromOptions(step);

    if(restart==0)
    {
        PetscViewerHDF5Open(PETSC_COMM_WORLD,fname,FILE_MODE_APPEND,&viewer);
        PetscViewerHDF5PushGroup(viewer,"/u_t");     

        VecCopy(u_0,u[0]);
        sprintf(dsname, "%08d",0);
        PetscObjectSetName((PetscObject)(u[0]),dsname);
        VecView(u[0],viewer); 
        while(its<maxIts)
        {
            MatMultAdd(A,u[its],b,u[its+1]);
            its++;
            sprintf(dsname, "%08d",its);
            PetscObjectSetName((PetscObject)(u[its]),dsname);
            VecView(u[its],viewer); 
        }
        VecSet(step,its);
        PetscObjectSetName((PetscObject)(step),"step");
        VecView(step,viewer); 
        PetscViewerHDF5PopGroup(viewer);
    }
    else
    {
        file_id=H5Fopen(fname,H5F_ACC_RDONLY,H5P_DEFAULT);
        group_id=H5Gopen(file_id,"/u_t",H5P_DEFAULT);
        dataset_id=H5Dopen(group_id,"step",H5P_DEFAULT);
        H5Dread(dataset_id,H5T_NATIVE_INT32,H5S_ALL,H5S_ALL,H5P_DEFAULT,&istep);
        H5Gclose(group_id);
        H5Fclose(file_id);      

        its=istep;
        sprintf(dsname, "%08d",its);
        PetscViewerHDF5Open(PETSC_COMM_WORLD,fname,FILE_MODE_UPDATE,&viewer);
        PetscViewerHDF5PushGroup(viewer,"/u_t");     
        PetscObjectSetName((PetscObject)(u[its]),dsname);
        VecLoad(u[its],viewer);

        while(its<maxIts)
        {
            MatMultAdd(A,u[its],b,u[its+1]);
            its++;
            sprintf(dsname, "%08d",its);
            PetscObjectSetName((PetscObject)(u[its]),dsname);
            VecView(u[its],viewer); 

        }
        VecSet(step,its);
        PetscObjectSetName((PetscObject)(step),"step");
        VecView(step,viewer); 

        PetscViewerHDF5PopGroup(viewer);
    }
    

    // MatView(A,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(u_0,PETSC_VIEWER_STDOUT_WORLD);

    VecDestroy(&u_0);
    PetscViewerDestroy(&viewer);
   
    PetscFinalize();
    return 0;

}
