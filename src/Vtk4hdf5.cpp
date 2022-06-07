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

#include "H5pubconf.h"
#include "vtk_hdf5.h"

int main(int argc, char * argv[])
{
  const std::string inputFN("u_t_output.hdf5");
  const std::string outputFN("u_t_visualization");
  const int         canreadfile;

  canreadfile = CanReadFile(inputFN.c_str());
  if(canreadfile = 1)
  {
    printf("file exits and can be read.")
  }
  else
  {
    printf("file can not be read.")
  }

  vtkHDFReader * reader = vtkHDFReader::New();

  reader -> Delete();

  printf("Here.")

  return 0;
}


// EOF
