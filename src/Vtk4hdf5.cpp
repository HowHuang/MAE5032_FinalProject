#include "vtkFloatArray.h"
#include "vtkPoints.h"

#include "vtk_hdf5.h"
#include "H5public.h"
#include <math.h>
#include <stdlib.h>

#define FILE "../default.hdf5"
#define GROUP "u_t"

int main(int argc, char * argv[])
{
  htri_t              file_test;
  hid_t               file_id, group_id, dataset_id;
  hsize_t             ds_num, ds_size;
  herr_t              status;
  char                dsname[10000];
  
  if((file_test = H5Fis_hdf5(FILE)) > 0)
  {
    printf("file type of temperature (u) for vtk visualization is committed\n");
  }
  else if (!file_test)
  {
    printf("file type of temperature (u) for vtk visualization should be hdf5\n");
    return -1;
  } 
  else
  {
    printf("error determining whether data type is committed\n");
    return -1;
  }

  file_id = H5Fopen(FILE,H5F_ACC_RDONLY,H5P_DEFAULT);
  /*
    H5F_ACC_RDWR:   Allow read and write access to file
    H5F_ACC_RDONLY: Allow read-only access to file
    */

  group_id = H5Gopen(file_id,GROUP,H5P_DEFAULT);
  if(group_id < 0)
  {
    printf("group type of temperature (u) in hdf5 file for vtk is invalid\n");
    return -1;
  }

  printf("group was opened.\n");

  status = H5Gget_num_objs(group_id,&ds_num);
  printf("The number of datasets in the group (u_t): %d\n", ds_num);
  /*
    another way to get the group's size (the number of dataset with the group):
    */
    // H5G_info_t          ginfo;
    // status = H5Gget_info(group_id, &ginfo);
    // printf("The number of datasets in the group (u_t): %d.\n", ginfo.nlinks);

  //use the first dataset to get the size of every dataset
  sprintf(dsname, "%08d", 0);
  dataset_id = H5Dopen(group_id,dsname,H5P_DEFAULT);
  ds_size = H5Dget_storage_size(dataset_id);
  printf("The number of data in each dataset: %d\n", ds_size/sizeof(double));
  /*
    another way to get the u_t's size:
    */
    //int       da_size;
    //da_size = ds_size/sizeof(double);
    //double    u_t[da_size];
    //std::cout << da_size << std::endl;
  
  double    *u_t = (double*)malloc(ds_size);
  // status = H5Dclose(dataset_id);               ///暂时注释，指读一个dataset测试

  //declare vkt objects
  int        n = sqrt(ds_size/sizeof(double));
  double     delta = 1.0/n;

  printf("The partitioning n for data is: %d*%d\n", n, n);

  vtkFloatArray * zvalues = vtkFloatArray::New();
  zvalues -> SetName("UValues");
  
  vtkPoints * points = vtkPoints::New();
  points -> SetDataTypeToDouble();

  for(int i=0; i<n; ++i)
  {
    for(int j=0; j<n; ++j)
    {
      points  -> InsertNextPoint(delta/2+i*delta,delta/2+i*delta,0);
      zvalues -> InsertNextValue(u_t[i*n+j]);
      // std::cout << delta/2+i*delta << std::endl;
    }
  }
  
  
  status = H5Dclose(dataset_id);


  // for(int img = 1; img < ds_num-1; ++img)
  // {
  //   sprintf(dsname, "%08d", img);
  //   dataset_id = H5Dopen(group_id,dsname,H5P_DEFAULT);

  //   ds_size = H5Dget_storage_size(dataset_id);
    
  //   status = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,u_t);



  //   status = H5Dclose(dataset_id);
  //   printf("dataset was opened.\n");
  // }
  


  

  points -> Delete();

  // 


  // H5Gget_info(group_id, &storage_type);
  // printf("group type is: %c.\n", * storage_type);

  // // group_num = H5Fget_obj_count(file_id, H5F_OBJ_FILE);
  // // printf("The number of group: %d", group_num);


  status = H5Gclose(group_id);
  status = H5Fclose(file_id);
  return 0;
}


// EOF
