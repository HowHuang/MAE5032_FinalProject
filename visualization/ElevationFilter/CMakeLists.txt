PROJECT(vtkHDF-leFilter)

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# SET(VTK_DIR /work/ese-huangh/lib/VTK-8.2.0/lib64/cmake/vtk-8.2)
SET(VTK_DIR /usr/local/lib/cmake/vtk-8.2)

# FIND_PACKAGE(VTK COMPONENTS vtkCommonCore vtkCommonSystem
#   vtkCommonMisc vtkCommonMath vtkIOCore vtkIOLegacy
#   vtkIOXML vtkIOParallelXML vtkhdf5 REQUIRED)
find_package(VTK COMPONENTS 
  vtkhdf5
  vtkCommonColor 
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonTransforms
  vtkFiltersCore
  vtkFiltersGeneral
  vtkFiltersGeometry
  vtkFiltersModeling
  vtkFiltersSources
  vtkIOCore vtkIOLegacy vtkIOXML
  vtkIOMovie
  vtkIOExportOpenGL2
  vtkIOExportPDF
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  REQUIRED
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "Visulisation: Unable to find the VTK build folder.\n")
  message(FATAL_ERROR "This software tested through VTK 8.2.0.")
endif()

INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS}/vtkhdf5)
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS}/vtkhdf5/src)
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS}/vtkhdf5/hl/src)
include(${VTK_USE_FILE})

ADD_EXECUTABLE( vtk_EleFilter Vtk4hdf5_ElevationFilter_MultiFiles.cpp )
TARGET_LINK_LIBRARIES( vtk_EleFilter ${VTK_LIBRARIES} )

# EOF