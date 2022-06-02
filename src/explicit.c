static char help[] = "Calculating a_{l,m}^{t+1} by the explicit method.";

#include <petsc.h>

static PetscErrorCode VecLoadIfExists_Private(Vec b,PetscViewer viewer,PetscBool *has)
{
  PetscBool      hdf5=PETSC_FALSE;

  PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERHDF5,&hdf5);
  if (hdf5) {
#if defined(PETSC_HAVE_HDF5)
  PetscViewerHDF5HasObject(viewer,(PetscObject)b,has);
  if (*has) VecLoad(b,viewer);
#else
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"PETSc must be configured with HDF5 to use this feature");
#endif
  } else {
    PetscErrorCode ierrp;
    PetscPushErrorHandler(PetscReturnErrorHandler,NULL);
    ierrp = VecLoad(b,viewer);
    PetscPopErrorHandler();
    *has  = ierrp ? PETSC_FALSE : PETSC_TRUE;
  }
  return 0;
}


int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank,size;
  PetscRandom    rand;
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
  PetscViewer    viewer;
  char           file[PETSC_MAX_PATH_LEN]="";    // input file name
  char           check_groupName;
  PetscBool      hdf5=PETSC_FALSE;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr); //size需要吗？
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  /*
    Determine reading the boundary data from the boundata.h5
    or randomly setting the boundary data according the dimension n.
    */
  ierr = PetscOptionsGetString(NULL,NULL,"-f",file,sizeof(file),NULL);CHKERRQ(ierr);
  /*
    Decide whether to use the HDF5 reader.
    */
  ierr = PetscOptionsGetBool(NULL,NULL,"-hdf5",&hdf5,NULL);CHKERRQ(ierr);
  
  /*
    Begin a segment of code that may be preloaded (run twice) to get accurate timings
    */
  ierr = PetscPreLoadBegin(PETSC_FALSE,"Load system");CHKERRQ(ierr); //需要预加载吗?
  /*
    Open hdf file.  Note that we use FILE_MODE_READ to indicate
    reading from this file. This code refers to:
    https://petsc.org/release/src/ksp/ksp/tutorials/ex27.c.html
    */
  if (hdf5) {
#if defined(PETSC_HAVE_HDF5)
    PetscViewerHDF5Open(PETSC_COMM_WORLD,file,FILE_MODE_READ,&viewer);
    PetscViewerPushFormat(viewer,PETSC_VIEWER_HDF5_MAT);
#else
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"PETSc must be configured with HDF5 to use this feature");
#endif
  } else {
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&viewer);
  }
  /*
    Load the vectors if it is present in the file, otherwise set up random data.
    */
  ierr = VecLoadIfExists_Private(g_b,viewer,&has);CHKERRQ(ierr);
  if (!has) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Failed to load boundary data, so use random data.\n");CHKERRQ(ierr);
    // randomly set boundary data，and the length is dependent on n
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rand);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rand);CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &g_b);CHKERRQ(ierr);
    ierr = VecSetSizes(x1, PETSC_DECIDE, n);
    ierr = VecSetFromOptions(g_b);CHKERRQ(ierr);
    ierr = VecSetRandom(g_b, rand);
    ierr = VecDuplicate(g_b, &g_t);CHKERRQ(ierr);
    ierr = VecDuplicate(g_b, &g_l);CHKERRQ(ierr);
    ierr = VecDuplicate(g_b, &g_r);CHKERRQ(ierr);
    ierr = VecDuplicate(g_b, &h_b);CHKERRQ(ierr);
    ierr = VecDuplicate(g_b, &h_t);CHKERRQ(ierr);
    ierr = VecDuplicate(g_b, &h_l);CHKERRQ(ierr);
    ierr = VecDuplicate(g_b, &h_r);CHKERRQ(ierr);
    //怎么确保g和h不同时为0，也不同时不为零

    //print the randomly set data
    PetscPrintf(PETSC_COMM_WORLD,"The randomly set boundary g data along bottom, top, left, right face'.\n");
    ierr = VecView(g_b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(g_t,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(g_l,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(g_r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"The randomly set boundary h data along bottom, top, left, right face'.\n");
    ierr = VecView(h_b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(h_t,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(h_l,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecView(h_r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  } else if {
    ierr = VecCreate(PETSC_COMM_WORLD,&g_b);CHKERRQ(ierr);
    ierr = VecSetFromOptions(g_b);CHKERRQ(ierr);
    //bottom and check the group name
    ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_u/bottom_u");
    ierr = PetscViewerHDF5GetGroup(viewer, check_groupName);
    if (check_groupName != "/boundary_u/bottom_u") {
      PetscPrintf(PETSC_COMM_WORLD,"Please use unified group name, which should be '/boundary_u/bottom_u'.\n");
      PetscPrintf(PETSC_COMM_WORLD,"'/boundary_u/top_u','/boundary_u/left_u','/boundary_u/right_u'\n");
    }
    ierr = VecLoad(g_b, viewer); //需要PetscObjectSetName((PetscObject) x2r, "x2");设名字吗？
    ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
    //top
    ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_u/top_u");
    ierr = VecDuplicate(g_b, &g_t);CHKERRQ(ierr);
    ierr = VecLoad(g_t, viewer);
    ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
    //left
    ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_u/left_u");
    ierr = VecDuplicate(g_b, &g_l);CHKERRQ(ierr);
    ierr = VecLoad(g_l, viewer);
    ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
    //right
    ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_u/left_u");
    ierr = VecDuplicate(g_b, &g_r);CHKERRQ(ierr);
    ierr = VecLoad(g_r, viewer);
    ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
    //similarly for h
    ierr = VecDuplicate(g_b, &h_b);CHKERRQ(ierr);
    //also check the group name
    ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_h/bottom_h");
    ierr = PetscViewerHDF5GetGroup(viewer, check_groupName);
    if (check_groupName != "/boundary_h/bottom_h") {
      PetscPrintf(PETSC_COMM_WORLD,"Please use unified group name, which should be '/boundary_h/bottom_h'.\n");
      PetscPrintf(PETSC_COMM_WORLD,"'/boundary_h/top_h','/boundary_h/left_h','/boundary_h/right_h'\n");
    }
    ierr = VecLoad(h_b, viewer); //需要PetscObjectSetName((PetscObject) x2r, "x2");设名字吗？
    ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
    ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_h/top_h");
    ierr = VecDuplicate(g_b, &h_t);CHKERRQ(ierr);
    ierr = VecLoad(h_t, viewer);
    ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
    ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_h/left_h");
    ierr = VecDuplicate(g_b, &h_l);CHKERRQ(ierr);
    ierr = VecLoad(h_l, viewer);
    ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
    ierr = PetscViewerHDF5PushGroup(viewer, "/boundary_h/left_h");
    ierr = VecDuplicate(g_b, &h_r);CHKERRQ(ierr);
    ierr = VecLoad(h_r, viewer);
    ierr = PetscViewerHDF5PopGroup(viewer);CHKERRQ(ierr);
  }
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);


  // Create u
  nn = (n+1)*(n+1); //size of u

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