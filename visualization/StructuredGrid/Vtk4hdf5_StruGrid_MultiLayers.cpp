#include "vtkDoubleArray.h"
#include "vtkArrayData.h"
#include "vtkCellData.h"
#include "vtkCellIterator.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkLookupTable.h"
#include "vtkTypedDataArray.h"
#include "vtkPolyDataMapper.h"
// #include "vtkRenderWindow.h"
// #include "vtkRenderWindowInteractor.h"
// #include "vtkRenderer.h"
#include "vtkNamedColors.h"
// #include "vtkCamera.h"
// #include "vtkWarpScalar.h"
// #include "vtkXMLPolyDataWriter.h"
// #include "vtkType.h"
#include "vtkStructuredGrid.h"
#include "vtkXMLStructuredGridWriter.h"

#include <iterator>
#include <map>
#include <set>

#include "vtk_hdf5.h"
#include "H5public.h"
#include <math.h>
#include <stdlib.h>

#define FILE "../../default.hdf5"
#define GROUP "u_t"

int main(int argc, char * argv[])
{
  htri_t              file_test;
  hid_t               file_id, group_id, dataset_id;
  hsize_t             ds_num, ds_size;
  herr_t              status;
  const std::string   filename("../output/StruGrid");
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
  printf("The number of datasets in the group (u_t): %lld\n", ds_num);
  /*
    another way to get the group's size (the number of dataset with the group):
    */
    // H5G_info_t          ginfo;
    // status = H5Gget_info(group_id, &ginfo);
    // printf("The number of datasets in the group (u_t): %d.\n", ginfo.nlinks);

  // use the first dataset to get the size of every dataset
  sprintf(dsname, "%08d", 0);
  dataset_id = H5Dopen(group_id,dsname,H5P_DEFAULT);
  ds_size    = H5Dget_storage_size(dataset_id);
  printf("The number of data in each dataset: %lld\n", ds_size/sizeof(double));
  /*
    another way to get the u_t's size:
    */
    //int       da_size;
    //da_size = ds_size/sizeof(double);
    //double    u_t[da_size];
    //std::cout << da_size << std::endl;
  
  double   *u_t = (double*)malloc(ds_size);
  status = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,u_t);
  status = H5Dclose(dataset_id);               ///暂时注释，指读一个dataset测试

  /*****************************
   *Start declaring vkt objects*
   *****************************/
  unsigned int      n = sqrt(ds_size/sizeof(double));
  double            delta = 1.0/n;
  auto              dataSize = n*n; //*(ds_num-1)
  auto              numberOfCells = (n-1)*(n-1);

  printf("The partitioning n for data is: %d*%d\n", n, n);

  for(int img = 0; img < ds_num-1; ++img)
  {
    sprintf(dsname, "%08d", img);
    dataset_id = H5Dopen(group_id,dsname,H5P_DEFAULT);
    ds_size    = H5Dget_storage_size(dataset_id);
    status     = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,u_t);
  
    vtkDoubleArray * pointValues = vtkDoubleArray::New();
    pointValues -> SetName("UValues");
    pointValues -> SetNumberOfComponents(1);
    pointValues -> SetNumberOfTuples(dataSize);
    for (size_t k=0; k<dataSize; ++k)
    {
      pointValues->SetValue(k, u_t[k]);
      // std::cout << u_t[k] << std::endl;
    }

    vtkDoubleArray * cellValues = vtkDoubleArray::New();
    cellValues  -> SetNumberOfTuples(numberOfCells);
    for (size_t k=0; k<numberOfCells; ++k)
    {
      cellValues->SetValue(k, u_t[k]);
    }

    vtkPoints * points = vtkPoints::New();

    for(size_t i=0; i<n; ++i)
    {
      for(size_t j=0; j<n; ++j)
      {
        points  -> InsertNextPoint(delta/2+i*delta,delta/2+j*delta,0);//0.1*img
        //if 3-D, the z values can be u_t[i*n+j]
      }
    }
    status = H5Dclose(dataset_id);
  
    // Create a grid
    vtkStructuredGrid * structuredGrid = vtkStructuredGrid::New();
    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(static_cast<int>(n), static_cast<int>(n),
                                  static_cast<int>(1));
    structuredGrid->SetPoints(points);
    structuredGrid->GetCellData()->SetScalars(cellValues);
    structuredGrid->GetPointData()->SetScalars(pointValues);

    // The key is the cell Id, the value is a set of corresponding point Ids.
    std::map<vtkIdType, std::set<vtkIdType>> cellPointIds;
    vtkCellIterator* it = structuredGrid->NewCellIterator();
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
    {
      vtkIdList * pointIds = it->GetPointIds();
      std::set<vtkIdType> ptIds;

      // std::cout << pointIds->GetId(0) << std::endl;
      vtkIdType  numId  = pointIds->GetNumberOfIds();
      // std::cout << "2::" << numId << std::endl;
      for (vtkIdType id = pointIds->GetId(0);id != pointIds->GetId(0)+pointIds->GetNumberOfIds();++id)
      {
        ptIds.insert(id);
      }
      cellPointIds[it->GetCellId()] = ptIds;
    }
    it->Delete();

    // std::cout << "Cells and their points" << std::endl;
    // for (auto cell : cellPointIds)
    // {
    //   std::cout << "Cell Id: " << cell.first << " Point Ids: ";
    //   for (auto id = cell.second.begin(); id != cell.second.end(); ++id)
    //     if (id != std::prev(cell.second.end()))
    //     {
    //       std::cout << *id << ", ";
    //     }
    //     else
    //     {
    //       std::cout << *id << std::endl;
    //     }
    // }

    // The key is the point Id and the value is a set of corresponding cell Ids.
    std::map<vtkIdType, std::set<vtkIdType>> commonPointIds;
    for (auto cell : cellPointIds)
    {
       for (auto pointId : cell.second)
      {
        commonPointIds[pointId].insert(cell.first);
      }
    }

    // std::cout << "Point Ids shared between cells" << std::endl;
    // for (auto point = commonPointIds.begin(); point != commonPointIds.end(); ++point)
    // {
    //   if (point->second.size() > 1)
    //   {
    //     std::cout << "Point Id: " << point->first << " CellIds: ";
    //     for (auto cellId = point->second.begin(); cellId != point->second.end(); ++cellId)
    //     {
    //       if (cellId != std::prev(point->second.end()))
    //       {
    //         std::cout << *cellId << ", ";
    //       }
    //       else
    //       {
    //         std::cout << *cellId << std::endl;
    //       }
    //     }
    //   }
    // }

  /*****************
   * Setup outputs *
   * ***************/
    vtkXMLStructuredGridWriter * writer = vtkXMLStructuredGridWriter::New();
    std::string name_to_write(filename);
    
    name_to_write.append("_t");
    name_to_write.append(std::to_string(img));
    name_to_write.append(".vts");
    writer -> SetFileName(name_to_write.c_str());
    writer -> SetInputData(structuredGrid);
    writer -> Write();

    printf("dataset %d was writen.\n", img);

    pointValues     -> Delete();
    points          -> Delete();
    cellValues      -> Delete();
    structuredGrid  -> Delete();
    writer          -> Delete();
  }
  
  status = H5Gclose(group_id);
  status = H5Fclose(file_id);
  std::cout << "Success!" << std::endl;
  return EXIT_SUCCESS;
}


// EOF
