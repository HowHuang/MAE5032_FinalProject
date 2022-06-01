static char help[] = "Calculating a_{l,m}^{t+1} by the explicit method.";

#include <petsc.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  PetscInt       i,j,istart,iend,n=10,nlocal;
  PetscInt       nn;
  PetscInt       col[3];

  PetscScalar    u0;
  Vec            u_t,u_tplus
  nn = 






  // Initialization
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  ierr = VecSetSizes(u,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(u,&istart,&iend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(u,&nlocal);CHKERRQ(ierr);

  if( rank == 0 )
  {
   for (i=1; i<n; i++)
    {
     for (j=1; j<n; j++) //u_{1 to 9, 1 to 9}
     {
      ierr = VecSetValues(u,1,&i,&u0,INSERT_VALUES);CHKERRQ(ierr);
     }
    }
  }

  





}