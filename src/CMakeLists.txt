PROJECT(vtk-HDF)

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

SET(HDF5_DIR /work/ese-huangh/lib/hdf5-1.12.1/lib)
SET(HDF5_LIBRARIES /work/ese-huangh/lib/hdf5-1.12.1/lib)
SET(HDF5_INCLUDE_DIRS /work/ese-huangh/lib/hdf5-1.12.1/include)
find_package(HDF5 COMPONENTS REQUIRED)

SET(VTK_DIR /work/ese-huangh/lib/VTK-8.2.0/lib64/cmake/vtk-8.2)
FIND_PACKAGE(VTK COMPONENTS vtkCommonCore vtkCommonSystem
  vtkCommonMisc vtkCommonMath vtkIOCore vtkIOLegacy
  vtkIOXML vtkIOParallelXML vtkIOHDF REQUIRED)

INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})

ADD_EXECUTABLE( vtkHDF Vtk4hdf5.cpp )
TARGET_LINK_LIBRARIES( vtkHDF ${HDF5_LIBRARIES} )
TARGET_LINK_LIBRARIES( vtkHDF ${VTK_LIBRARIES} )

# EOF