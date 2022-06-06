#include <petsc.h>
#include <hdf5.h>

typedef struct InputPara
{
    PetscScalar dt;     PetscScalar dl;     PetscScalar rho;
    PetscScalar c;      PetscScalar k;      PetscScalar f;
}InputPara;

inline PetscScalar CalElementOfRightVec(InputPara *IP, PetscScalar U1, PetscInt U2, PetscScalar H1, PetscScalar H2)
{   
    // PetscScalar P1 = IP->f * IP->dt / (IP->rho * IP->c);
    // PetscScalar P2_1 = U1 * 2 * IP->k * (pow(IP->dt,2)) / (pow(IP->rho,2) * pow(IP->c,2) * (pow(IP->dl,4)));
    // PetscScalar P2_2 = U2 * 2 * IP->k * (pow(IP->dt,2)) / (pow(IP->rho,2) * pow(IP->c,2) * (pow(IP->dl,4)));
    // PetscScalar P3_1 = H1 * IP->dt / (IP->rho * IP->c * pow(IP->dl,2));
    // PetscScalar P3_2 = H2 * IP->dt / (IP->rho * IP->c * pow(IP->dl,2));

    return  IP->f * IP->dt / (IP->rho * IP->c) + \
            U1 * 2 * IP->k * (pow(IP->dt,2)) / (pow(IP->rho,2) * pow(IP->c,2) * (pow(IP->dl,4))) + \
            U2 * 2 * IP->k * (pow(IP->dt,2)) / (pow(IP->rho,2) * pow(IP->c,2) * (pow(IP->dl,4))) + \
            H1 * IP->dt / (IP->rho * IP->c * pow(IP->dl,2)) + \
            H2 * IP->dt / (IP->rho * IP->c * pow(IP->dl,2));
}
    
int main(int argc,char **argv)
{
    //Vec             g_b,g_t,g_l,g_r;    //DIM = n x 1
    //Vec             h_b,h_t,h_l,h_r;    //DIM = n x 1
    Vec             u_0,b,u_t,u_tplus,temp;  //DIM = (n*n) x 1
    Mat             A;                  //DIM = (n*n) x (n*n)
    PetscViewer     viewerI,viewerO;    //I for input, O for output
    PetscMPIInt     rank,size;
    PetscInt        i, j, r, n = 10;
    MPI_Comm        comm;
    PetscErrorCode  ierr;
    PetscScalar     bb;
    PetscChar       ifname[PETSC_MAX_PATH_LEN]="g_fixed.hdf5";
    PetscChar       ofname[PETSC_MAX_PATH_LEN]="u_t_output.hdf5";
    PetscChar       dsname[PETSC_MAX_PATH_LEN]="default";
    PetscInt        col[5];
    PetscScalar     value[5];
    PetscScalar     W,N,E,S,P;
    PetscScalar     U1, U2, H1, H2;
    InputPara       IP; //PetscScalar     dt,dl,rho,c,k,f;

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

    
    W=(IP.k*IP.dt)/(IP.rho*IP.c*PetscPowScalar(IP.dl,2));
    N=(IP.k*IP.dt)/(IP.rho*IP.c*PetscPowScalar(IP.dl,2));
    E=(IP.k*IP.dt)/(IP.rho*IP.c*PetscPowScalar(IP.dl,2));
    S=(IP.k*IP.dt)/(IP.rho*IP.c*PetscPowScalar(IP.dl,2));
    P=1-(4*IP.k*IP.dt)/(IP.rho*IP.c*PetscPowScalar(IP.dl,2));
    // PetscPrintf(comm,"dt:%g, dl: %g, rho:%g, c:%g, k:%d\n",IP.dt,IP.dl,IP.rho,IP.k,IP.c);
    // PetscPrintf(comm,"W:%g,N:%g,E:%g,S:%g,P:%g\n",W,N,E,S,P);

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
    /*
    // PetscViewerHDF5PushGroup(viewerI,"/boundary");

    // VecCreate(comm,&g_b);
    // PetscObjectSetName((PetscObject)g_b, "g_bottom");
    // VecCreate(comm,&h_b);
    // PetscObjectSetName((PetscObject)h_b, "h_bottom");
    
    // VecCreate(comm,&g_t);
    // PetscObjectSetName((PetscObject)g_t, "g_top");
    // VecCreate(comm,&h_t);
    // PetscObjectSetName((PetscObject)h_t, "h_top");
    
    // VecCreate(comm,&g_r);
    // PetscObjectSetName((PetscObject)g_r, "g_right");  
    // VecCreate(comm,&h_r);
    // PetscObjectSetName((PetscObject)h_r, "h_right");
    
    // VecCreate(comm,&g_l);
    // PetscObjectSetName((PetscObject)g_l, "g_left"); 
    // VecCreate(comm,&h_l);
    // PetscObjectSetName((PetscObject)h_l, "h_left");

    // VecLoad(g_b,viewerI);
    // VecLoad(h_b,viewerI);

    // VecLoad(g_t,viewerI);
    // VecLoad(h_t,viewerI);

    // VecLoad(g_r,viewerI);
    // VecLoad(h_r,viewerI);

    // VecLoad(g_l,viewerI);
    // VecLoad(h_l,viewerI);        

    // PetscViewerHDF5PopGroup(viewerI);
    */
    
    PetscViewerHDF5PushGroup(viewerI,"/u0_2D");
    VecCreate(comm,&u_0);
    PetscObjectSetName((PetscObject)u_0, "u0_2DInit");
    VecLoad(u_0,viewerI);
    PetscViewerHDF5PopGroup(viewerI);

    PetscViewerDestroy(&viewerI);
    // ~ Load finished

    // //~ Print vecs
    // VecView(g_b,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(h_b,PETSC_VIEWER_STDOUT_WORLD);     
    // VecView(g_t,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(h_t,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(g_r,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(h_r,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(g_l,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(h_l,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(u_0,PETSC_VIEWER_STDOUT_WORLD);  

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
                        col[0]=r;col[1]=r+1;col[2]=r+n;
                        value[0]=P;value[1]=E;value[2]=S;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);


                        U1=g_t[j];  // VecGetValues(g_t,1,&j,&U1);
                        H1=h_t[j];  // VecGetValues(h_t,1,&j,&H1);
                        U2=g_l[i];  // VecGetValues(g_l,1,&i,&U2);
                        H2=h_l[i];  // VecGetValues(h_l,1,&i,&H2);
                        bb = CalElementOfRightVec(&IP,U1,U2,H1,H2);
                        VecSetValue(b,r,bb,INSERT_VALUES);
                        
                    }
                    else if(j==n-1)     //! 右上角的点
                    {
                        col[0]=r-1;col[1]=r;col[2]=r+n;
                        value[0]=W;value[1]=P;value[2]=S;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);

                        U1=g_t[j];  // VecGetValues(g_t,1,&j,&U1);
                        H1=h_t[j];  // VecGetValues(h_t,1,&j,&H1);
                        U2=g_r[i];  // VecGetValues(g_r,1,&i,&U2);
                        H2=h_r[i];  // VecGetValues(h_r,1,&i,&H2);             
                        bb = CalElementOfRightVec(&IP,U1,U2,H1,H2);
                        VecSetValue(b,r,bb,INSERT_VALUES);       
                    }
                    else                //! 上边界的点（除顶点外）
                    {
                        col[0]=r-1;col[1]=r;col[2]=r+1;col[3]=r+n;
                        value[0]=W;value[1]=P;value[2]=E;value[3]=S;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);

                        U1=g_t[j];  // VecGetValues(g_t,1,&j,&U1);
                        H1=h_t[j];  // VecGetValues(h_t,1,&j,&H1);
                        U2=0.0;
                        H2=0.0;
                        bb = CalElementOfRightVec(&IP,U1,U2,H1,H2);
                        VecSetValue(b,r,bb,INSERT_VALUES);              
                    }
                }
                else if(i==n-1)
                {
                    if(j==0)            //! 左下角的点
                    {
                        col[0]=r-n;col[1]=r;col[2]=r+1;
                        value[0]=N;value[1]=P;value[2]=E;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);  

                        U1=g_b[j];  // VecGetValues(g_b,1,&j,&U1);
                        H1=h_b[j];  // VecGetValues(h_b,1,&j,&H1);
                        U2=g_l[i];  // VecGetValues(g_l,1,&i,&U2);
                        H2=h_l[i];  // VecGetValues(h_l,1,&i,&H2);
                        bb = CalElementOfRightVec(&IP,U1,U2,H1,H2);
                        VecSetValue(b,r,bb,INSERT_VALUES);                                   
                    }
                    else if(j==n-1)     //! 右下角的点
                    {
                        col[0]=r-n;col[1]=r-1;col[2]=r;
                        value[0]=N;value[1]=W;value[2]=P;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);      

                        U1=g_b[i];  // VecGetValues(g_b,1,&j,&U1);
                        H1=h_b[j];  // VecGetValues(h_b,1,&j,&H1);
                        U2=g_r[i];  // VecGetValues(g_r,1,&i,&U2);
                        H2=h_r[i];  // VecGetValues(h_r,1,&i,&H2);
                        bb = CalElementOfRightVec(&IP,U1,U2,H1,H2);
                        VecSetValue(b,r,bb,INSERT_VALUES);                                     
                    }
                    else                //! 下边界的点 （除顶点外）
                    {
                        col[0]=r-n;col[1]=r-1;col[2]=r;col[3]=r+1;
                        value[0]=N;value[1]=W;value[2]=P;value[3]=E;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);  

                        U1=g_b[j];  // VecGetValues(g_b,1,&j,&U1);
                        H1=h_b[j];  // VecGetValues(h_b,1,&j,&H1);
                        U2=0.0;
                        H2=0.0;
                        bb = CalElementOfRightVec(&IP,U1,U2,H1,H2);
                        VecSetValue(b,r,bb,INSERT_VALUES);                                       
                    }
                }
                else
                {
                    if(j==0)            //! 左边界的点（除顶点外）
                    {
                        col[0]=r-n;col[1]=r;col[2]=r+1;col[3]=r+n;
                        value[0]=N;value[1]=P;value[2]=E;value[3]=S;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);
                        col[0]=r-n;col[1]=r;col[2]=r+1;
                        value[0]=N;value[1]=P;value[2]=E;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);  

                        U1=g_l[i];  // VecGetValues(g_l,1,&i,&U1);
                        H1=h_l[i];// VecGetValues(h_l,1,&i,&H1);
                        U2=0.0;
                        H2=0.0;
                        bb = CalElementOfRightVec(&IP,U1,U2,H1,H2);
                        VecSetValue(b,r,bb,INSERT_VALUES);

                    }
                    else if(j==n-1)     //! 右边界的点（除顶点外）
                    {
                        col[0]=r-n;col[1]=r-1;col[2]=r;col[3]=r+n;
                        value[0]=N;value[1]=W;value[2]=P;value[3]=S;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);

                        U1=g_r[i];  // VecGetValues(g_r,1,&i,&U1);
                        H1=h_r[i];  // VecGetValues(h_r,1,&i,&H1);
                        U2=0.0;
                        H2=0.0;
                        bb = CalElementOfRightVec(&IP,U1,U2,H1,H2);
                        VecSetValue(b,r,bb,INSERT_VALUES);                  
                    }
                    else                //! 内部点
                    {
                        col[0]=r-n;col[1]=r-1;col[2]=r;col[3]=r+1;col[4]=r+n;
                        value[0]=N;value[1]=W;value[2]=P;value[3]=E;value[4]=S;
                        MatSetValues(A,1,&r,5,col,value,INSERT_VALUES);

                        U1=0.0;
                        H1=0.0;
                        U2=0.0;
                        H2=0.0;
                        bb = CalElementOfRightVec(&IP,U1,U2,H1,H2);
                        VecSetValue(b,r,bb,INSERT_VALUES);
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
