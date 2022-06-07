PROJECT(vtk-HDF)

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

SET(VTK_DIR /work/ese-huangh/lib/VTK-8.2.0/lib64/cmake/vtk-8.2)
FIND_PACKAGE(VTK COMPONENTS vtkCommonCore vtkCommonSystem
  vtkCommonMisc vtkCommonMath vtkIOCore vtkIOLegacy
  vtkIOXML vtkIOParallelXML vtkhdf5 REQUIRED)

INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})

ADD_EXECUTABLE( vtkHDF Vtk4hdf5.cpp )
TARGET_LINK_LIBRARIES( vtkHDF ${VTK_LIBRARIES} )

# EOF