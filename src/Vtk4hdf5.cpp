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

int main(int argc, char * argv[])
{
  htri_t      file_test;
  hid_t       file_id, group_id, dataset_id;  //identifiers

  herr_t      group_name, status;

  
  if ((file_test = H5Fis_hdf5(FILE)) > 0)
  {
    printf("data type is committed\n");
  }
  else if (!file_test)
  {
    printf("data type is not committed\n");
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

  printf("file opened.");

  group_name = 

  //group_id = H5Gopen(file_id,"/boundary",H5P_DEFAULT);

  status = H5Fclose(file_id); 
  return 0;
}


// EOF
