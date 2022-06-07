#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkTetra.h"
#include "vtkGenericCell.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLGenericDataObjectReader.h"

#include "vtk_hdf5.h"
#include "H5public.h"

#define FILE "u_t_output.hdf5"
#define GROUP "u_t"

int main(int argc, char * argv[])
{
  htri_t       file_test;
  hid_t        file_id, group_id, dataset_id;  //identifiers
  hid_t        group_info;
  ssize_t      group_num;
  herr_t       group_status, status;
  H5G_info_t * storage_type;

  
  if((file_test = H5Fis_hdf5(FILE)) > 0)
  {
    printf("file type of temperature (u) for vtk visualization is committed\n");
  }
  else if (!file_test)
  {
    printf("file type of temperature (u) for vtk visualization should be hdf5\n");
  } 
  else
  {
    printf("error determining whether data type is committed\n");
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
  }

  printf("group was opened.");


  

  // H5Gget_info(group_id, &storage_type);
  // printf("group type is: %c.\n", * storage_type);

  // // group_num = H5Fget_obj_count(file_id, H5F_OBJ_FILE);
  // // printf("The number of group: %d", group_num);


  status = H5Fclose(file_id); 
  return 0;
}


// EOF
