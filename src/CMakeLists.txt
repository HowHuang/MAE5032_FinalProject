PROJECT(vtk-demo)

CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

SET(MPI_DIR /share/intel/2020u4/compilers_and_libraries_2020.4.304/linux/mpi)
find_package(MPI REQUIRED)
SET(VTK_DIR /work/ese-huangh/lib/VTK-8.2.0/lib64/cmake/vtk-8.2)
FIND_PACKAGE(VTK COMPONENTS vtkCommonCore vtkCommonSystem
  vtkCommonMisc vtkCommonMath vtkIOCore vtkIOLegacy vtkIOXML vtkIOParallelXML REQUIRED)

INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})

ADD_EXECUTABLE( vtkTESTmpi t6_cube4tet_TESTmpi.cpp )
TARGET_LINK_LIBRARIES( vtkTESTmpi ${VTK_LIBRARIES} )
target_link_libraries( vtkTESTmpi ${MPI_LIBRARIES} )

# EOF
