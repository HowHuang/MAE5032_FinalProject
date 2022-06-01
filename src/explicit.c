static char help[] = "Calculating a_{l,m}^{t+1} by the explicit method.";

#include <petsc.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank,size;
  PetscInt       i,j,istart,iend,n=10,nlocal;
  PetscInt       nn;
  PetscInt       col[3];

  PetscScalar    u0,kappa,rho,c;

  /*
    u_t refers to data including internal and boundary nodes.
    */
  Vec            u_t,u_tplus;
  /*
    _b: the bottom face; _t: the top face;
    _l: the left face;   _r: the right face;
    These four vecs will load from the corresponding .h5 files.
    in the data folder.
    */
  Vec            g_b,g_t,g_l,g_r;
  Vec            h_b,h_t,h_l,h_r;

  /*
    Method of reading the vecs refers to:
    https://petsc.org/release/src/vec/vec/tutorials/ex19.c.html
    We may be able to add a code for testing I/O, like:
    https://petsc.org/release/src/vec/vec/tutorials/ex10.c.html 
    */
  PetscViewer    fd;
  char           file[PETSC_MAX_PATH_LEN]="";    // input file name
  PetscBool      hdf5=PETSC_FALSE;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  /*
    Determine reading the boundary data from the boundata.h5
    or randomly setting the boundary data according the dimension n.
    */
  ierr = PetscOptionsGetString(NULL,NULL,"-f",file,sizeof(file),NULL);CHKERRQ(ierr);
  /*
    Decide whether to use the HDF5 reader.
    */
  ierr = PetscOptionsGetBool(NULL,NULL,"-hdf5",&hdf5,NULL);CHKERRQ(ierr);

  ierr = PetscPreLoadBegin(PETSC_FALSE,"Load system");CHKERRQ(ierr);
  /*
    Open hdf file.  Note that we use FILE_MODE_READ to indicate
    reading from this file. This code refers to:
    https://petsc.org/release/src/ksp/ksp/tutorials/ex27.c.html
    */
  if (hdf5) {
#if defined(PETSC_HAVE_HDF5)
    PetscViewerHDF5Open(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);
    PetscViewerPushFormat(fd,PETSC_VIEWER_HDF5_MAT);
#else
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"PETSc must be configured with HDF5 to use this feature");
#endif
  } else {
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);
  }


  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "boundata.h5", FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&g_b);CHKERRQ(ierr);
  ierr = VecSetSizes(g_b,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(g_b);CHKERRQ(ierr);//Still need this?




  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);


  nn = (n+1)*(n+1); //size of u

  // Loading all boundary vectors






  // Initialization of u
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  ierr = VecSetSizes(u,PETSC_DECIDE,nn);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(u,&istart,&iend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(u,&nlocal);CHKERRQ(ierr);

  if( rank == 0 )
  {
   // Internal nodes
   for (i=1; i<n; i++)
    {
     for (j=1; j<n; j++) //u_{1 to 9, 1 to 9}
     {
      ierr = VecSetValues(u,1,&i,&u0,INSERT_VALUES);CHKERRQ(ierr);
     }
    }
   // Boundary nodes
   i = 0;

  }

  // 
  




}