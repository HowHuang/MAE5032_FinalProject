static char help[] = "Calculating a_{l,m}^{t+1} by the explicit method.";

#include <petsc.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  PetscInt       i,j,istart,iend,n=10,nlocal;
  PetscInt       nn;
  PetscInt       col[3];

  PetscScalar    u0,kappa,rho,c;

  /*
    u_t refers to data including internal and boundary nodes.
    */
  Vec            u_t,u_tplus;
  /*
    u_b: the bottom face; u_t: the top face;
    u_l: the left face; u_r: the right face;
    These four vecs will load from the corresponding .h5 files.
    in the data folder.
    */
  Vec            u_b,u_t,u_l,u_r;

  /*
    method of reading the vecs refers to:
    https://petsc.org/release/src/vec/vec/tutorials/ex10.c.html
    */
  PetscViewer       viewer;
  PetscBool         vstage2,vstage3,mpiio_use,isbinary = PETSC_FALSE;
#if defined(PETSC_HAVE_HDF5)
  PetscBool         ishdf5 = PETSC_FALSE;
#endif
#if defined(PETSC_HAVE_ADIOS)
  PetscBool         isadios = PETSC_FALSE;
#endif
  PetscScalar const *values;
#if defined(PETSC_USE_LOG)
  PetscLogEvent  VECTOR_GENERATE,VECTOR_READ;
#endif

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  mpiio_use = vstage2 = vstage3 = PETSC_FALSE;
  PetscOptionsGetBool(NULL,NULL,"-binary",&isbinary,NULL);
#if defined(PETSC_HAVE_HDF5)
  PetscOptionsGetBool(NULL,NULL,"-hdf5",&ishdf5,NULL);
#endif
#if defined(PETSC_HAVE_ADIOS)
  PetscOptionsGetBool(NULL,NULL,"-adios",&isadios,NULL);
#endif
  PetscOptionsGetBool(NULL,NULL,"-mpiio",&mpiio_use,NULL);
  PetscOptionsGetBool(NULL,NULL,"-sizes_set",&vstage2,NULL);
  PetscOptionsGetBool(NULL,NULL,"-type_set",&vstage3,NULL);

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
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