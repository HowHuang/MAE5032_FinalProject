#include <petsc.h>
#include <hdf5.h>
#include "explicit.h"
/*
typedef struct InputPara
{
    PetscScalar dt;     PetscScalar dl;     PetscScalar rho;
    PetscScalar c;      PetscScalar k;      PetscScalar f;
}InputPara;

typedef struct Bound
{
    PetscScalar ub;     PetscScalar hb;
    PetscScalar ut;     PetscScalar ht;
    PetscScalar ul;     PetscScalar hl;
    PetscScalar ur;     PetscScalar hr;    
}Bound; // ~ for a point

typedef struct IterMaterial
{
    PetscScalar W;      PetscScalar E;
    PetscScalar N;      PetscScalar S;
    PetscScalar P;      PetscScalar b;
}IterMaterial;

enum Location
{
    RightTop    = 1,    TopSide     = 2,
    LeftTop     = 3,    LeftSide    = 4,
    LeftBottom  = 5,    BottomSide  = 6,
    RightBottom = 7,    RightSide   = 8,
    Internal    = 9
};
*/
void CalIterationMaterial(InputPara* IP, Bound* bound, IterMaterial* IM, enum Location loc)
{
    //~ InputPara       Known;
    //~ bound           Known;
    //~ IterMaterial    Unknown;
    PetscScalar base = (IP->k * IP->dt) / (IP->rho * IP->c * PetscPowScalar(IP->dl,2));
    PetscScalar partF = IP->f * IP->dt / (IP->rho * IP->c);
    PetscScalar partH = bound->hb * IP->dt / (IP->rho * IP->c * pow(IP->dl,2)) + \
                        bound->ht * IP->dt / (IP->rho * IP->c * pow(IP->dl,2)) + \
                        bound->hl * IP->dt / (IP->rho * IP->c * pow(IP->dl,2)) + \
                        bound->hr * IP->dt / (IP->rho * IP->c * pow(IP->dl,2));
    PetscScalar partU = 2 * base * (bound->ut + bound->ur + bound->ub + bound->ul);                    
    switch (loc)
    {
    case RightTop:
    {
        IM->W = base;
        IM->S = base;
        IM->b = partF + partU + partH;            
        if(bound->ut!=0 && bound->ur!=0 && bound->ht==0 && bound->hr==0)
            IM->P = 1-(1+1+2+2)*base;
        else if(bound->ut==0 && bound->ur!=0 && bound->ht==0 && bound->hr==0)
            IM->P = 1-(1+1+2)*base;
        else if(bound->ut==0 && bound->ur!=0 && bound->ht==0 && bound->hr==0)
            IM->P = 1-(1+1+2)*base;
        else if(bound->ut==0 && bound->ur==0 && bound->ht==0 && bound->hr==0)
            IM->P = 1-(1+1+2+2)*base;
        else if (bound->ut==0 && bound->ur==0 && (bound->ht!=0 || bound->hr==0))
            IM->P = 1-(1+1)*base;
        else printf("u and h can not be both given at same boundary");
    }
    break;

    case LeftTop:
    {
        IM->E = base;
        IM->S = base;
        IM->b = partF + partU + partH;            
        if(bound->ut!=0 && bound->ul!=0 && bound->ht==0 && bound->hl==0)
            IM->P = 1-(1+1+2+2)*base;
        else if(bound->ut==0 && bound->ul!=0 && bound->ht==0 && bound->hl==0)
            IM->P = 1-(1+1+2)*base;
        else if(bound->ut==0 && bound->ul!=0 && bound->ht==0 && bound->hl==0)
            IM->P = 1-(1+1+2)*base;
        else if(bound->ut==0 && bound->ul==0 && bound->ht==0 && bound->hl==0)
            IM->P = 1-(1+1+2+2)*base;
        else if (bound->ut==0 && bound->ul==0 && (bound->ht!=0 || bound->hl==0))
            IM->P = 1-(1+1)*base;
        else printf("u and h can not be both given at same boundary");
    }
    break;

    case LeftBottom:
    {
        IM->E = base;
        IM->S = base;
        IM->b = partF + partU + partH;            
        if(bound->ub!=0 && bound->ul!=0 && bound->hb==0 && bound->hl==0)
            IM->P = 1-(1+1+2+2)*base;
        else if(bound->ub==0 && bound->ul!=0 && bound->hb==0 && bound->hl==0)
            IM->P = 1-(1+1+2)*base;
        else if(bound->ub==0 && bound->ul!=0 && bound->hb==0 && bound->hl==0)
            IM->P = 1-(1+1+2)*base;
        else if(bound->ub==0 && bound->ul==0 && bound->hb==0 && bound->hl==0)
            IM->P = 1-(1+1+2+2)*base;
        else if (bound->ub==0 && bound->ul==0 && (bound->hb!=0 || bound->hl==0))
            IM->P = 1-(1+1)*base;
        else printf("u and h can not be both given at same boundary");
    }
    break;

    case RightBottom:
    {
        IM->E = base;
        IM->S = base;
        IM->b = partF + partU + partH;            
        if(bound->ub!=0 && bound->ur!=0 && bound->hb==0 && bound->hr==0)
            IM->P = 1-(1+1+2+2)*base;
        else if(bound->ub==0 && bound->ur!=0 && bound->hb==0 && bound->hr==0)
            IM->P = 1-(1+1+2)*base;
        else if(bound->ub==0 && bound->ur!=0 && bound->hb==0 && bound->hr==0)
            IM->P = 1-(1+1+2)*base;
        else if(bound->ub==0 && bound->ur==0 && bound->hb==0 && bound->hr==0)
            IM->P = 1-(1+1+2+2)*base;
        else if (bound->ub==0 && bound->ur==0 && (bound->hb!=0 || bound->hr==0))
            IM->P = 1-(1+1)*base;
        else printf("u and h can not be both given at same boundary");
    }
    break;

    case RightSide:
    {
        IM->W = base;
        IM->S = base;
        IM->N = base;
        IM->b = partF + partU + partH;            
        if(bound->hr==0)
            IM->P = 1-(1+1+1+2)*base;
        else if (bound->hr!=0 && bound->ur==0)
            IM->P = 1-(1+1+1)*base;
        else printf("u and h can not be both given at same boundary\n");
    }
    break;

    case TopSide:
    {
        IM->W = base;
        IM->E = base;
        IM->S = base;
        IM->b = partF + partU + partH;   
        if(bound->ht==0)
            IM->P = 1-(1+1+1+2)*base;
        else if(bound->ht!=0 && bound->ut==0)
            IM->P = 1-(1+1+1)*base;
        else printf("u and h can not be both given at same boundary\n");
    }
    break;

    case LeftSide:
    {
        IM->N = base;
        IM->E = base;
        IM->S = base;
        IM->b = partF + partU + partH;   
        if(bound->hl==0)
            IM->P = 1-(1+1+1+2)*base;
        else if(bound->hl!=0 && bound->ul==0)
            IM->P = 1-(1+1+1)*base;
        else printf("u and h can not be both given at same boundary\n");
    }   
    break;

    case BottomSide:
    {
        IM->N = base;
        IM->E = base;
        IM->W = base;
        IM->b = partF + partU + partH;   
        if(bound->hb==0)
            IM->P = 1-(1+1+1+2)*base;
        else if(bound->hb!=0 && bound->ub==0)
            IM->P = 1-(1+1+1)*base;
        else printf("u and h can not be both given at same boundary\n");
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
    Vec             u_0,b,u_t,u_tplus,temp;  //DIM = (n*n) x 1
    Mat             A;                  //DIM = (n*n) x (n*n)
    PetscViewer     viewerI,viewerO;    //I for input, O for output
    PetscMPIInt     rank,size;
    PetscInt        i, j, r, n = 10;
    MPI_Comm        comm;
    PetscErrorCode  ierr;
    PetscChar       ifname[PETSC_MAX_PATH_LEN]="g_fixed.hdf5";
    PetscChar       ofname[PETSC_MAX_PATH_LEN]="u_t_output.hdf5";
    PetscChar       dsname[PETSC_MAX_PATH_LEN]="default";
    PetscInt        col[5];
    PetscScalar     value[5];

    InputPara       IP; //PetscScalar     dt,dl,rho,c,k,f;
    Bound           bound;
    IterMaterial    IM;
    Location        loc;

    PetscScalar     *g_b, *g_t, *g_l, *g_r;
    PetscScalar     *h_b, *h_t, *h_l, *h_r;

    // iteration parts
    PetscInt        its=0 , maxIts = 100;
    Vec             *u;         

    ierr = PetscInitialize(&argc, &argv, (char*) 0, NULL);if (ierr) return ierr;
    comm = PETSC_COMM_WORLD;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
    ierr = MPI_Comm_size(comm, &size);CHKERRMPI(ierr);

    ierr = PetscOptionsGetInt(NULL,NULL, "-n", &n, NULL);
    ierr = PetscOptionsGetString(NULL,NULL,"-ifname",ifname,sizeof(ifname),NULL);
    ierr = PetscOptionsGetString(NULL,NULL,"-ofname",ofname,sizeof(ofname),NULL);
    ierr = PetscOptionsGetInt(NULL,NULL, "-maxIts", &maxIts, NULL);

    ierr = PetscOptionsGetScalar(NULL,NULL,"-dt",  &(IP.dt),  NULL);
    ierr = PetscOptionsGetScalar(NULL,NULL,"-dl",  &(IP.dl),  NULL);
    ierr = PetscOptionsGetScalar(NULL,NULL,"-rho", &(IP.rho), NULL);
    ierr = PetscOptionsGetScalar(NULL,NULL,"-c",   &(IP.c),   NULL); 
    ierr = PetscOptionsGetScalar(NULL,NULL,"-k",   &(IP.k),   NULL);
    ierr = PetscOptionsGetScalar(NULL,NULL,"-f",   &(IP.f),   NULL);

    PetscMalloc4(n, &g_b, n, &g_t, n, &g_l, n, &g_r);
    PetscMalloc4(n, &h_b, n, &h_t, n, &h_l, n, &h_r);

    hid_t   file_id, group_id, dataset_id;
    file_id=H5Fopen(ifname,H5F_ACC_RDONLY,H5P_DEFAULT);
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

    H5Gclose(group_id);
    H5Fclose(file_id);

    // ~ Load data from specified hdf5 file.
    PetscViewerHDF5Open(PETSC_COMM_WORLD,ifname,FILE_MODE_READ,&viewerI);
    
    PetscViewerHDF5PushGroup(viewerI,"/u0_2D");
    VecCreate(comm,&u_0);
    PetscObjectSetName((PetscObject)u_0, "u0_2DInit");
    VecLoad(u_0,viewerI);
    PetscViewerHDF5PopGroup(viewerI);

    PetscViewerDestroy(&viewerI);
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
                        bound.ut=g_t[j];bound.ht=h_t[j];
                        bound.ul=g_l[i];bound.hl=h_l[i];
                        bound.ub=0;     bound.hb=0;
                        bound.ur=0;     bound.hr=0;
                        CalIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r;col[1]=r+1;col[2]=r+n;
                        value[0]=IM.P;value[1]=IM.E;value[2]=IM.S;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES); 
                        VecSetValue(b,r,IM.b,INSERT_VALUES);
                    }
                    else if(j==n-1)     //! 右上角的点
                    {
                        loc = RightTop;
                        bound.ut=g_t[j];bound.ht=h_t[j];
                        bound.ur=g_r[i];bound.hr=h_r[i];
                        bound.ub=0;     bound.hb=0;
                        bound.ul=0;     bound.hl=0;
                        CalIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-1;col[1]=r;col[2]=r+n;
                        value[0]=IM.W;value[1]=IM.P;value[2]=IM.S;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);   
                        VecSetValue(b,r,IM.b,INSERT_VALUES);       
                    }
                    else                //! 上边界的点（除顶点外）
                    {   
                        loc = TopSide;
                        bound.ut=g_t[j];bound.ht=h_t[j];
                        bound.ur=0;     bound.hr=0;
                        bound.ub=0;     bound.hb=0;
                        bound.ul=0;     bound.hl=0;
                        CalIterationMaterial(&IP,&bound,&IM,loc);
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
                        bound.ub=g_b[j];bound.hb=h_b[j];
                        bound.ur=0;     bound.hr=0;
                        bound.ut=0;     bound.ht=0;
                        bound.ul=g_l[i];bound.hl=h_l[i];
                        CalIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-n;col[1]=r;col[2]=r+1;
                        value[0]=IM.N;value[1]=IM.P;value[2]=IM.E;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);  
                        VecSetValue(b,r,IM.b,INSERT_VALUES);                                   
                    }
                    else if(j==n-1)     //! 右下角的点
                    {   
                        loc = RightBottom;
                        bound.ub=g_b[j];bound.hb=h_b[j];
                        bound.ul=0;     bound.hl=0;
                        bound.ut=0;     bound.ht=0;
                        bound.ur=g_r[i];bound.hr=h_r[i];
                        CalIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-n;col[1]=r-1;col[2]=r;
                        value[0]=IM.N;value[1]=IM.W;value[2]=IM.P;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);      
                        VecSetValue(b,r,IM.b,INSERT_VALUES);                                     
                    }
                    else                //! 下边界的点 （除顶点外）
                    {
                        loc = BottomSide;
                        bound.ub=g_b[j];bound.hb=h_b[j];
                        bound.ul=0;     bound.hl=0;
                        bound.ut=0;     bound.ht=0;
                        bound.ur=0;     bound.hr=0;
                        CalIterationMaterial(&IP,&bound,&IM,loc);
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
                        bound.ub=0;     bound.hb=0;
                        bound.ul=g_l[i];bound.hl=h_l[i];
                        bound.ut=0;     bound.ht=0;
                        bound.ur=0;     bound.hr=0;
                        CalIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-n;col[1]=r;col[2]=r+1;col[3]=r+n;
                        value[0]=IM.N;value[1]=IM.P;value[2]=IM.E;value[3]=IM.S;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);
                        VecSetValue(b,r,IM.b,INSERT_VALUES);

                    }
                    else if(j==n-1)     //! 右边界的点（除顶点外）
                    {
                        loc = RightSide;
                        bound.ub=0;     bound.hb=0;
                        bound.ur=g_r[i];bound.hr=h_r[i];
                        bound.ut=0;     bound.ht=0;
                        bound.ul=0;     bound.hl=0;
                        CalIterationMaterial(&IP,&bound,&IM,loc);
                        col[0]=r-n;col[1]=r-1;col[2]=r;col[3]=r+n;
                        value[0]=IM.N;value[1]=IM.W;value[2]=IM.P;value[3]=IM.S;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);
                        VecSetValue(b,r,IM.b,INSERT_VALUES);                  
                    }
                    else                //! 内部点
                    {
                        loc = Internal;
                        bound.ub=0;     bound.hb=0;
                        bound.ur=0;     bound.hr=0;
                        bound.ut=0;     bound.ht=0;
                        bound.ul=0;     bound.hl=0;
                        CalIterationMaterial(&IP,&bound,&IM,loc);
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

    PetscViewerHDF5Open(PETSC_COMM_WORLD,ofname,FILE_MODE_WRITE,&viewerO);
    PetscViewerHDF5PushGroup(viewerO,"/u_t");
    
    VecCopy(u_0,u[0]);
    sprintf(dsname, "%d",0);
    PetscObjectSetName((PetscObject)(u[0]),dsname);
    VecView(u[0],viewerO); 
    while(its<maxIts)
    {
        MatMultAdd(A,u[its],b,u[its+1]);
        its++;
        sprintf(dsname, "%d",its);
        PetscObjectSetName((PetscObject)(u[its]),dsname);
        VecView(u[its],viewerO); 
    }

    PetscViewerHDF5PopGroup(viewerO);

    MatView(A,PETSC_VIEWER_STDOUT_WORLD);
    VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    VecView(u_0,PETSC_VIEWER_STDOUT_WORLD);

    VecDestroy(&u_0);
    PetscViewerDestroy(&viewerI);
   
    PetscFinalize();
    return 0;

}
