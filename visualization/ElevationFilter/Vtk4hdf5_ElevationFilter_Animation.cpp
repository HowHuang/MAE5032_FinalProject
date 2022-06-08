#include "vtkFloatArray.h"
#include "vtkArrayData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkDelaunay2D.h"
#include "vtkElevationFilter.h"
#include "vtkLookupTable.h"
#include "vtkTypedDataArray.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkNamedColors.h"
#include "vtkCamera.h"
#include "vtkWarpScalar.h"
// #include "vtkXMLPolyDataWriter.h"
#include "vtkActor.h"
#include "vtkWindowToImageFilter.h"
#include "vtkGenericMovieWriter.h"
#include "vtkSmartPointer.h"

#include "vtk_hdf5.h"
#include "H5public.h"
#include <math.h>
#include <stdlib.h>

#define FILE "../../default.hdf5"
#define GROUP "u_t"

/*
  https://gitlab.kitware.com/vtk/vtk/-/tree/v8.2.0/IO/Movie
  why my 8.2 does not have "vtkAVIWriter.h"????
  another try: vtkOggTheoraWriter.h??
  */

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
  sprintf(dsname, "%08d", 1);
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
  status = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,u_t);
  // status = H5Dclose(dataset_id);               ///暂时注释，指读一个dataset测试

  //declare vkt objects
  int        n = sqrt(ds_size/sizeof(double));
  double     delta = 1.0/n;

  printf("The partitioning n for data is: %d*%d\n", n, n);

  // vtkFloatArray * zvalues = vtkFloatArray::New();
  // zvalues -> SetName("UValues");
  
  vtkPoints * points = vtkPoints::New();
  points -> SetDataTypeToDouble();

  for(int i=0; i<n; ++i)
  {
    for(int j=0; j<n; ++j)
    {
      points  -> InsertNextPoint(delta/2+i*delta,delta/2+j*delta,u_t[i*n+j]);
      // zvalues -> InsertNextValue(u_t[i*n+j]);
      // std::cout << u_t[i*n+j] << std::endl;
    }
  }
  status = H5Dclose(dataset_id);
  
  double    bounds[6];
  points -> GetBounds(bounds);
  /*GetBounds() Return a pointer to the geometry bounding box in the form 
    (xmin,xmax, ymin, ymax, zmin, zmax).
    */
  printf("The min and max of u : %f and %f\n", bounds[4], bounds[5]);

  // Add the grid points to a polydata object
  vtkPolyData * inputPolyData = vtkPolyData::New();
  inputPolyData -> SetPoints(points);

  // Triangulate the grid points
  vtkDelaunay2D * delaunay = vtkDelaunay2D::New();
  delaunay -> SetInputData(inputPolyData);
  delaunay -> Update();

  vtkElevationFilter * elevationFilter = vtkElevationFilter::New();
  elevationFilter -> SetInputConnection(delaunay->GetOutputPort());
  elevationFilter -> SetLowPoint(0.0, 0.0, bounds[4]);
  elevationFilter -> SetHighPoint(0.0, 0.0, bounds[5]);
  elevationFilter -> Update();

  vtkPolyData * output = vtkPolyData::New();
  output -> ShallowCopy(dynamic_cast<vtkPolyData*>(elevationFilter->GetOutput()));
  std::cout << "output here: " << output->GetNumberOfPoints() << std::endl;

  vtkFloatArray * elevation = vtkFloatArray::New();
  elevation = dynamic_cast<vtkFloatArray*>(output->GetPointData()->GetArray("Elevation"));
    /*GetArray will Returns the n-th vtkArray in the collection.
      https://vtk.org/doc/nightly/html/classvtkArrayData.html#ae404ca06bbffe9927ac59cb908c4ab95
      */

  // Create the color map
  vtkLookupTable * colorLookupTable = vtkLookupTable::New();
  colorLookupTable -> SetTableRange(bounds[4], bounds[5]);
  colorLookupTable -> Build();

  // Generate the colors for each point based on the color map
  vtkUnsignedCharArray * colors = vtkUnsignedCharArray::New();
  colors -> SetNumberOfComponents(3);
  colors -> SetName("Colors");

  for (vtkIdType i = 0; i < output->GetNumberOfPoints(); i++)
  {
    double val = elevation->GetValue(i);
    // std::cout << "val: " << val << std::endl; //checked

    double dcolor[3];
    colorLookupTable -> GetColor(val, dcolor);
    // std::cout << "dcolor: " << dcolor[0] << " " << dcolor[1] << " " <<
    // dcolor[2] << std::endl;
    unsigned char color[3];
    for (unsigned int j = 0; j < 3; j++)
    {
      color[j] = 255 * dcolor[j] / 1.0;
    }
    // std::cout << "color: " << (int)color[0] << " " << (int)color[1] << " " <<
    // (int)color[2] << std::endl;

    colors -> InsertNextTypedTuple(color);
    /*
      InsertNextTupleValue is now InsertNextTypedTuple:
      https://vtk.org/doc/nightly/html/VTK-7-1-Changes.html
      */
  }

  output -> GetPointData() -> AddArray(colors);
  colors -> Delete();

  /***********************
   * Setup visualization *
   * *********************/
  vtkNamedColors * Namecolors = vtkNamedColors::New();
  // Set the background color.
  vtkColor3d backgroundColor = Namecolors->GetColor3d("SlateGray");
  std::cout << "here" << std::endl;
  std::cout << output << std::endl;

  // Create a mapper and actor
  vtkPolyDataMapper * mapper = vtkPolyDataMapper::New();
  vtkActor          * actor  = vtkActor::New();
  mapper -> SetInputData(output);
  actor  -> SetMapper(mapper);
  // actor  -> GetProperty()->SetColor(Namecolors->GetColor3d("Cyan").GetData());
  std::cout << "here3" << std::endl;

  // Create a renderer, render window, and interactor
  vtkRenderer     * renderer     = vtkRenderer::New();
  vtkRenderWindow * renderWindow = vtkRenderWindow::New();
  renderWindow -> AddRenderer(renderer);
  // renderWindow -> SetSize(640, 480);
  renderWindow -> SetWindowName("DataAnimation");
/*
  vtkRenderWindowInteractor * renderWindowInteractor = vtkRenderWindowInteractor::New();
  renderWindowInteractor -> SetRenderWindow(renderWindow);

  // Initialize must be called prior to creating timer events.
  renderWindowInteractor -> Initialize();
  renderWindowInteractor -> CreateRepeatingTimer(ds_num-1);
*/

  // vtkCallbackCommand * timerCallback = vtkCallbackCommand::New();
  // timerCallback -> SetCallback(TimerCallbackFunction);
  // timerCallback->SetClientData(programmableFilter);

  // Add the actor to the scene
  renderer -> AddActor(actor);
  renderer -> SetBackground(backgroundColor.GetData());

  // Render and interact
  renderWindow->Render();
  auto camera = renderer->GetActiveCamera();
  camera->SetPosition(2.26841, -1.51874, 1.805);
  camera->SetFocalPoint(-0.148582, 0.0814323, 0.24803);
  camera->SetViewUp(0.157813, 0.800687, 0.577923);
  camera->SetDistance(3.29037);
  camera->SetClippingRange(1.14823, 5.60288);
  // or
  // renderWindow->Render();
  // renderer->GetActiveCamera()->Azimuth(30);
  // renderer->GetActiveCamera()->Elevation(30);
  // renderer->ResetCameraClippingRange();



  /*****************
   * Setup outputs *
   * ***************/

  vtkGenericMovieWriter * writer;
    // vtkNew<vtkGenericMovieWriter> writer;
  // vtkSmartPointer<vtkAVIWriter> writer=vtkSmartPointer<vtkAVIWriter>::New();

  vtkWindowToImageFilter * filter = vtkWindowToImageFilter::New();
  filter -> SetInput(renderWindow);
  writer -> SetInputConnection(filter->GetOutputPort());
  writer -> SetFileName("demo_EleFilter.avi");
  // writer -> InlineDataOn();
  // writer -> SetRenderWindow(renderWindow);
  writer -> Write();

  writer -> Delete();

  // for(int img = 1; img < ds_num-1; ++img)
  // {
  //   sprintf(dsname, "%08d", img);
  //   dataset_id = H5Dopen(group_id,dsname,H5P_DEFAULT);

  //   ds_size = H5Dget_storage_size(dataset_id);
    
  //   status = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,u_t);



  //   status = H5Dclose(dataset_id);
  //   printf("dataset was opened.\n");
  // }
  


  
  // zvalues -> Delete();
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