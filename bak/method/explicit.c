static char help[] = "Calculating a_{l,m}^{t+1} by the explicit method.";

#include <petsc.h>

int main(int *argc,char ***argv)
{
  
  PetscErrorCode ierr;
  PetscMPIInt    rank,size;
  PetscRandom    rand;
  PetscInt       i,j,istart,iend,n=10,nlocal;
  PetscInt       nn;

  PetscScalar    u0,kappa,rho,c;
  Vec            u;
  Vec            u_t,u_tplus;     //temperature to be calculated
  
  Vec            g_b,g_t,g_l,g_r; //boundary g
  Vec            h_b,h_t,h_l,h_r; //boundray h

  // Method of reading the vecs refers to:
  // https://petsc.org/release/src/vec/vec/tutorials/ex19.c.html
  // We may be able to add a code for testing I/O, like:
  // https://petsc.org/release/src/vec/vec/tutorials/ex10.c.html 

  PetscViewer    viewer;
  PetscChar      file[PETSC_MAX_PATH_LEN]="";    // input file name
  PetscChar      check_groupName;
  PetscBool      has;

  //assert(!PETSC_HAVE_HDF5);
  #ifndef PETSC_HAVE_HDF5
  printf("Error: PETSc must be configured with HDF5 to use this feature.\n");
  return -1;
  #endif

  // initialize PETSc and MPI
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

  // Determine reading the boundary data from the boundata.h5
  // or randomly setting the boundary data according the dimension n.
  // 用户指定边界条件数据文件的地址
  ierr = PetscOptionsGetString(NULL,NULL,"-f",file,sizeof(file),NULL);CHKERRQ(ierr);
  
  // Open hdf file.  Note that we use FILE_MODE_READ to indicate
  // reading from this file. This code refers to:
  // https://petsc.org/release/src/ksp/ksp/tutorials/ex27.c.html

  PetscViewerHDF5Open(PETSC_COMM_WORLD,file,FILE_MODE_READ,&viewer);
  PetscViewerPushFormat(viewer,PETSC_VIEWER_HDF5_MAT);

  // Load the vectors if it is present in the file, otherwise set up random data.

  // ierr = VecLoadIfExists_Private(g_b,viewer,&has);CHKERRQ(ierr);
  // if (!has) 
  // {
  //   ierr = PetscPrintf(PETSC_COMM_WORLD,"Failed to load boundary data, so use random data.\n");CHKERRQ(ierr);
  //   // randomly set boundary data，and the length is dependent on n
  //   ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rand);CHKERRQ(ierr);
  //   ierr = PetscRandomSetFromOptions(rand);CHKERRQ(ierr);
  //   ierr = VecCreate(PETSC_COMM_WORLD, &g_b);CHKERRQ(ierr);
  //   ierr = VecSetSizes(g_b, PETSC_DECIDE, n);
  //   ierr = VecSetFromOptions(g_b);CHKERRQ(ierr);
  //   ierr = VecSetRandom(g_b, rand);
  //   ierr = VecDuplicate(g_b, &g_t);CHKERRQ(ierr);
  //   ierr = VecDuplicate(g_b, &g_l);CHKERRQ(ierr);
  //   ierr = VecDuplicate(g_b, &g_r);CHKERRQ(ierr);
  //   ierr = VecDuplicate(g_b, &h_b);CHKERRQ(ierr);
  //   ierr = VecDuplicate(g_b, &h_t);CHKERRQ(ierr);
  //   ierr = VecDuplicate(g_b, &h_l);CHKERRQ(ierr);
  //   ierr = VecDuplicate(g_b, &h_r);CHKERRQ(ierr);
  //   //怎么确保g和h不同时为0，也不同时不为零

  //   //print the randomly set data
  //   PetscPrintf(PETSC_COMM_WORLD,"The randomly set boundary g data along bottom, top, left, right face'.\n");
  //   ierr = VecView(g_b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //   ierr = VecView(g_t,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //   ierr = VecView(g_l,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //   ierr = VecView(g_r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //   PetscPrintf(PETSC_COMM_WORLD,"The randomly set boundary h data along bottom, top, left, right face'.\n");
  //   ierr = VecView(h_b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //   ierr = VecView(h_t,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //   ierr = VecView(h_l,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //   ierr = VecView(h_r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  // } 
  // else
  // {
  //   ierr = VecCreate(PETSC_COMM_WORLD,&g_b);CHKERRQ(ierr);
  //   ierr = VecSetFromOptions(g_b);CHKERRQ(ierr);
  //   //bottom and check the group name
  //   ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_u/bottom_u");
  //   ierr = PetscViewerHDF5GetGroup(viewer, check_groupName);
  //   if (check_groupName != "/boundary_u/bottom_u") {
  //     PetscPrintf(PETSC_COMM_WORLD,"Please use unified group name, which should be '/boundary_u/bottom_u'.\n");
  //     PetscPrintf(PETSC_COMM_WORLD,"'/boundary_u/top_u','/boundary_u/left_u','/boundary_u/right_u'\n");
  //   }
  //   ierr = VecLoad(g_b, viewer); //需要PetscObjectSetName((PetscObject) x2r, "x2");设名字吗？
  //   ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
  //   //top
  //   ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_u/top_u");
  //   ierr = VecDuplicate(g_b, &g_t);CHKERRQ(ierr);
  //   ierr = VecLoad(g_t, viewer);
  //   ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
  //   //left
  //   ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_u/left_u");
  //   ierr = VecDuplicate(g_b, &g_l);CHKERRQ(ierr);
  //   ierr = VecLoad(g_l, viewer);
  //   ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
  //   //right
  //   ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_u/left_u");
  //   ierr = VecDuplicate(g_b, &g_r);CHKERRQ(ierr);
  //   ierr = VecLoad(g_r, viewer);
  //   ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
  //   //similarly for h
  //   ierr = VecDuplicate(g_b, &h_b);CHKERRQ(ierr);
  //   //also check the group name
  //   ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_h/bottom_h");
  //   ierr = PetscViewerHDF5GetGroup(viewer, check_groupName);
  //   if (check_groupName != "/boundary_h/bottom_h") {
  //     PetscPrintf(PETSC_COMM_WORLD,"Please use unified group name, which should be '/boundary_h/bottom_h'.\n");
  //     PetscPrintf(PETSC_COMM_WORLD,"'/boundary_h/top_h','/boundary_h/left_h','/boundary_h/right_h'\n");
  //   }
  //   ierr = VecLoad(h_b, viewer); //需要PetscObjectSetName((PetscObject) x2r, "x2");设名字吗？
  //   ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
  //   ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_h/top_h");
  //   ierr = VecDuplicate(g_b, &h_t);CHKERRQ(ierr);
  //   ierr = VecLoad(h_t, viewer);
  //   ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
  //   ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_h/left_h");
  //   ierr = VecDuplicate(g_b, &h_l);CHKERRQ(ierr);
  //   ierr = VecLoad(h_l, viewer);
  //   ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
  //   ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_h/left_h");
  //   ierr = VecDuplicate(g_b, &h_r);CHKERRQ(ierr);
  //   ierr = VecLoad(h_r, viewer);
  //   ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
  // }
  // ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);


  // // Create u
  // nn = (n+1)*(n+1); //size of u

  // ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  // ierr = VecSetSizes(u,PETSC_DECIDE,nn);CHKERRQ(ierr);
  // ierr = VecSetFromOptions(u);CHKERRQ(ierr);

  // ierr = VecGetOwnershipRange(u,&istart,&iend);CHKERRQ(ierr);
  // ierr = VecGetLocalSize(u,&nlocal);CHKERRQ(ierr);

  // if( rank == 0 )
  // {
  //  // Internal nodes
  //  for (i=1; i<n; i++)
  //   {
  //    for (j=1; j<n; j++) //u_{1 to 9, 1 to 9}
  //    {
  //     ierr = VecSetValues(u,1,&i,&u0,INSERT_VALUES);CHKERRQ(ierr);
  //    }
  //   }
  //  // Boundary nodes
  //  i = 0;

  // }

  // // 
  

return 0;
}