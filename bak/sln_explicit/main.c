#include <petsc.h>

int main(int argc,char **argv)
{
    Vec             g_b,g_t,g_l,g_r;    //DIM = n x 1
    Vec             h_b,h_t,h_l,h_r;    //DIM = n x 1
    Vec             u_0,b,u_t,u_tplus;  //DIM = (n*n) x 1
    Mat             A;                  //DIM = (n*n) x (n*n)
    PetscInt        low,high,ldim,lsize;
    PetscInt        ilocal,iglobal;
    PetscViewer     viewerI,viewerO;    //I for input, O for output
    PetscMPIInt     rank,size;
    PetscInt        i, j, r, n = 10;
    MPI_Comm        comm;
    PetscErrorCode  ierr;
    PetscScalar     zero=0.0,one=1.0;
    PetscChar       ifname[PETSC_MAX_PATH_LEN]="g_fixed.hdf5";
    PetscChar       ofname[PETSC_MAX_PATH_LEN]="default_Out.hdf5";
    PetscInt        col[5];
    PetscScalar     value[5];
    PetscScalar     W,N,E,S,P;
    PetscScalar     dt,dl,rho,c,k;

    ierr = PetscInitialize(&argc, &argv, (char*) 0, NULL);if (ierr) return ierr;
    comm = PETSC_COMM_WORLD;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
    ierr = MPI_Comm_size(comm, &size);CHKERRMPI(ierr);

    ierr = PetscOptionsGetInt(NULL,NULL, "-n", &n, NULL);
    ierr = PetscOptionsGetString(NULL,NULL,"-ifname",ifname,sizeof(ifname),NULL);
    ierr = PetscOptionsGetString(NULL,NULL,"-ofname",ofname,sizeof(ofname),NULL);

    ierr = PetscOptionsGetScalar(NULL,NULL,"-dt",  &dt,  NULL);
    ierr = PetscOptionsGetScalar(NULL,NULL,"-dl",  &dl,  NULL);
    ierr = PetscOptionsGetScalar(NULL,NULL,"-rho", &rho, NULL);
    ierr = PetscOptionsGetScalar(NULL,NULL,"-c",   &c,   NULL); 
    ierr = PetscOptionsGetScalar(NULL,NULL,"-k",   &k,   NULL);
    
    W=(k*dt)/(rho*c*PetscPowScalar(dl,2));
    N=(k*dt)/(rho*c*PetscPowScalar(dl,2));
    E=(k*dt)/(rho*c*PetscPowScalar(dl,2));
    S=(k*dt)/(rho*c*PetscPowScalar(dl,2));
    P=1-(4*k*dt)/(rho*c*PetscPowScalar(dl,2));
    PetscPrintf(comm,"dt:%g, dl: %g, rho:%g, c:%g, k:%d\n",dt,dl,rho,k,c);
    PetscPrintf(comm,"W:%g,N:%g,E:%g,S:%g,P:%g\n",W,N,E,S,P);

    // ~ Load data from specified hdf5 file.
    PetscViewerHDF5Open(PETSC_COMM_WORLD,ifname,FILE_MODE_READ,&viewerI);
    PetscViewerHDF5PushGroup(viewerI,"/boundary");

    VecCreate(comm,&g_b);
    PetscObjectSetName((PetscObject)g_b, "g_bottom");
    VecCreate(comm,&h_b);
    PetscObjectSetName((PetscObject)h_b, "h_bottom");
    
    VecCreate(comm,&g_t);
    PetscObjectSetName((PetscObject)g_t, "g_top");
    VecCreate(comm,&h_t);
    PetscObjectSetName((PetscObject)h_t, "h_top");
    
    VecCreate(comm,&g_r);
    PetscObjectSetName((PetscObject)g_r, "g_right");  
    VecCreate(comm,&h_r);
    PetscObjectSetName((PetscObject)h_r, "h_right");
    
    VecCreate(comm,&g_l);
    PetscObjectSetName((PetscObject)g_l, "g_left"); 
    VecCreate(comm,&h_l);
    PetscObjectSetName((PetscObject)h_l, "h_left");

    VecLoad(g_b,viewerI);
    VecLoad(h_b,viewerI);

    VecLoad(g_t,viewerI);
    VecLoad(h_t,viewerI);

    VecLoad(g_r,viewerI);
    VecLoad(h_r,viewerI);

    VecLoad(g_l,viewerI);
    VecLoad(h_l,viewerI);

    PetscViewerHDF5PopGroup(viewerI);

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
    MatSetSizes(A,PETSC_DECIDE,n*n,n*n,n*n);
    MatSetFromOptions(A);
    MatSetUp(A);
    
    // ~ Generate the vector b for Au_t+b=u_tplus
    VecCreate(comm,&b);
    VecSetSizes(b,n,n*n);
    VecSetFromOptions(b);
    VecSetUp(b);

    // ~ Assign values to A and b

    if(rank==0)
    {
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                r=i*n+j;    //global row index of A
                if(i==0)
                {
                    if(j==0)
                    {
                        col[0]=r;col[1]=r+1;col[2]=r+n;
                        value[0]=P;value[1]=E;value[2]=S;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);
                    }
                    else if(j==n-1)
                    {
                        col[0]=r-1;col[1]=r;col[2]=r+n;
                        value[0]=W;value[1]=P;value[2]=S;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);
                    }
                    else
                    {
                        col[0]=r-1;col[1]=r;col[2]=r+1;col[3]=r+n;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);
                    }
                }
                else if(i==n-1)
                {
                    if(j==0)
                    {
                        col[0]=r-n;col[1]=r;col[2]=r+1;
                        value[0]=N;value[1]=P;value[2]=E;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);                        
                    }
                    else if(j==n-1)
                    {
                        col[0]=r-n;col[1]=r-1;col[2]=r;
                        value[0]=N;value[1]=W;value[2]=P;
                        MatSetValues(A,1,&r,3,col,value,INSERT_VALUES);                        
                    }
                    else
                    {
                        col[0]=r-n;col[1]=r-1;col[2]=r;col[3]=r+1;
                        value[0]=N;value[1]=W;value[2]=P;value[3]=E;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);                        
                    }
                }
                else
                {
                    if(j==0)
                    {
                        col[0]=r-n;col[1]=r;col[2]=r+1;col[3]=r+n;
                        value[0]=N;value[1]=P;value[2]=E;value[3]=S;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);
                    }
                    else if(j==n-1)
                    {
                        col[0]=r-n;col[1]=r-1;col[2]=r;col[3]=r+n;
                        value[0]=N;value[1]=W;value[2]=P;value[3]=S;
                        MatSetValues(A,1,&r,4,col,value,INSERT_VALUES);
                    }
                    else
                    {
                        col[0]=r-n;col[1]=r-1;col[2]=r;col[3]=r+1;col[4]=r+n;
                        value[0]=N;value[1]=W;value[2]=P;value[3]=E;value[4]=S;
                        MatSetValues(A,1,&r,5,col,value,INSERT_VALUES);
                    }
                }
            }
        }
    }

    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); 

    

    MatView(A,PETSC_VIEWER_STDOUT_WORLD);

    VecDestroy(&g_b);VecDestroy(&g_t);VecDestroy(&g_r);VecDestroy(&g_l);
    VecDestroy(&h_b);VecDestroy(&h_t);VecDestroy(&h_r);VecDestroy(&h_l);
    VecDestroy(&u_0);
    PetscViewerDestroy(&viewerI);
   
    PetscFinalize();
    return 0;

}
