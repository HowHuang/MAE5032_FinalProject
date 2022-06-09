PROJECT(vtkHDF-StruGrid)

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
  vtkFiltersCore
  vtkFiltersGeneral
  vtkFiltersGeometry
  vtkFiltersModeling
  vtkFiltersSources
  vtkIOCore vtkIOLegacy vtkIOXML vtkIOParallelXML
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  REQUIRED
)

INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS}/vtkhdf5)
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS}/vtkhdf5/src)
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS}/vtkhdf5/hl/src)

ADD_EXECUTABLE( vtk_StruGrid Vtk4hdf5_StruGrid_MultiLayers.cpp )
TARGET_LINK_LIBRARIES( vtk_StruGrid ${VTK_LIBRARIES} )

# EOF