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

void add_int_PointData( vtkPointSet * const &grid_w,
    const std::vector<int> &ptdata, const std::string &dataname );

void add_int_CellData( vtkPointSet * const &grid_w,
    const std::vector<int> &cldata, const std::string &dataname );

int main(int argc, char * argv[])
{
  const bool isXML = false;
  const std::string filename("cube2tet");
  int        n = 40;
  int        nn = (n+1)*(n+1)*(n+1);
  double     delta = 1.0/n;

  printf("point: %d\n", nn);
  printf("point: %f\n", delta);

  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  vtkPoints * ppt = vtkPoints::New();
  ppt->SetDataTypeToDouble();

  for(int i=0; i<(n+1); ++i)
  {
    for(int j=0; j<(n+1); ++j)
    {
      for(int k=0; k<(n+1); ++k)
      {
        ppt -> InsertPoint(i*(n+1)*(n+1)+j*(n+1)+k, 0.0+i*delta, 0.0+j*delta, 0.0+k*delta);
        printf("point: %-4d", i*(n+1)*(n+1)+j*(n+1)+k);
        printf("%f,%f,%f ",0.0+i*delta,0.0+j*delta,0.0+k*delta);
      }
    }
  }
  printf("here\n");
  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  vtkCell * cl = vtkTetra::New();
  for(int ii=0; ii<n; ++ii)
  {
    for(int jj=0; jj<n; ++jj)
    {
      for(int kk=0; kk<n; ++kk)
      {
        cl->GetPointIds()->SetId(0,  ii   *(n+1)*(n+1)+ jj   *(n+1)+kk  );
        cl->GetPointIds()->SetId(1,  ii   *(n+1)*(n+1)+ jj   *(n+1)+kk+1);
        cl->GetPointIds()->SetId(2, (ii+1)*(n+1)*(n+1)+ jj   *(n+1)+kk+1);
        cl->GetPointIds()->SetId(3,  ii   *(n+1)*(n+1)+(jj+1)*(n+1)+kk  );
        grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );
        printf("cube: %-4d", ii*n*n*6+jj*n*6+kk*6+0);

        cl->GetPointIds()->SetId(0,  ii   *(n+1)*(n+1)+ jj   *(n+1)+kk  );
        cl->GetPointIds()->SetId(1, (ii+1)*(n+1)*(n+1)+ jj   *(n+1)+kk  );
        cl->GetPointIds()->SetId(2, (ii+1)*(n+1)*(n+1)+ jj   *(n+1)+kk+1);
        cl->GetPointIds()->SetId(3,  ii   *(n+1)*(n+1)+(jj+1)*(n+1)+kk  );
        grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );
        printf("cube: %-4d", ii*n*n*6+jj*n*6+kk*6+1);

        cl->GetPointIds()->SetId(0, (ii+1)*(n+1)*(n+1)+ jj   *(n+1)+kk  );
        cl->GetPointIds()->SetId(1, (ii+1)*(n+1)*(n+1)+(jj+1)*(n+1)+kk  );
        cl->GetPointIds()->SetId(2, (ii+1)*(n+1)*(n+1)+ jj   *(n+1)+kk+1);
        cl->GetPointIds()->SetId(3,  ii   *(n+1)*(n+1)+(jj+1)*(n+1)+kk  );
        grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );
        printf("cube: %-4d", ii*n*n*6+jj*n*6+kk*6+2);

        cl->GetPointIds()->SetId(0, (ii+1)*(n+1)*(n+1)+(jj+1)*(n+1)+kk+1);
        cl->GetPointIds()->SetId(1, (ii+1)*(n+1)*(n+1)+(jj+1)*(n+1)+kk  );
        cl->GetPointIds()->SetId(2, (ii+1)*(n+1)*(n+1)+ jj   *(n+1)+kk+1);
        cl->GetPointIds()->SetId(3,  ii   *(n+1)*(n+1)+(jj+1)*(n+1)+kk  );
        grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );
        printf("cube: %-4d", ii*n*n*6+jj*n*6+kk*6+3);

        cl->GetPointIds()->SetId(0, (ii+1)*(n+1)*(n+1)+(jj+1)*(n+1)+kk+1);
        cl->GetPointIds()->SetId(1,  ii   *(n+1)*(n+1)+(jj+1)*(n+1)+kk+1);
        cl->GetPointIds()->SetId(2, (ii+1)*(n+1)*(n+1)+ jj   *(n+1)+kk+1);
        cl->GetPointIds()->SetId(3,  ii   *(n+1)*(n+1)+(jj+1)*(n+1)+kk  );
        grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );
        printf("cube: %-4d", ii*n*n*6+jj*n*6+kk*6+4);

        cl->GetPointIds()->SetId(0,  ii   *(n+1)*(n+1)+ jj   *(n+1)+kk+1);
        cl->GetPointIds()->SetId(1,  ii   *(n+1)*(n+1)+(jj+1)*(n+1)+kk+1);
        cl->GetPointIds()->SetId(2, (ii+1)*(n+1)*(n+1)+ jj   *(n+1)+kk+1);
        cl->GetPointIds()->SetId(3,  ii   *(n+1)*(n+1)+(jj+1)*(n+1)+kk  );
        grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );
        printf("cube: %-4d", ii*n*n*6+jj*n*6+kk*6+5);
      }
    }
  }

  cl -> Delete();
  printf("here2");

  std::vector<int> node_id = {};
  for(int nid=1; nid<(nn+1); ++nid)
  {
    node_id.push_back(nid);
  }
  add_int_PointData(grid_w, node_id, "GlobalNodeID");

  std::vector<int> cell_id = {};
  for(int cid=1; cid<(6*n*n*n+1); ++cid)
  {
    cell_id.push_back(cid);
  }
  add_int_CellData(grid_w, cell_id, "GlobalCellID");

  if( isXML )
  {
    vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
    std::string name_to_write(filename);
    name_to_write.append(".vtu");
    writer -> SetFileName( name_to_write.c_str() );

    writer->SetInputData(grid_w);
    writer->Write();
    writer->Delete();
  }
  else
  {
    vtkUnstructuredGridWriter * writer = vtkUnstructuredGridWriter::New();
    std::string name_to_write(filename);
    name_to_write.append(".vtk");
    writer -> SetFileName( name_to_write.c_str() );

    writer->SetInputData(grid_w);
    writer->Write();
    writer->Delete();
  }

  grid_w->Delete();
  return 0;
}

void add_int_PointData( vtkPointSet * const &grid_w,
    const std::vector<int> &ptdata, const std::string &dataname )
{
  vtkIntArray * data = vtkIntArray::New();
  data -> SetNumberOfComponents(1);
  data -> SetName(dataname.c_str());

  for(unsigned int ii=0; ii<ptdata.size(); ++ii)
    data -> InsertComponent(ii, 0, ptdata[ii]);

  grid_w -> GetPointData() -> AddArray( data );
  data -> Delete();
}

void add_int_CellData( vtkPointSet * const &grid_w,
    const std::vector<int> &cldata, const std::string &dataname )
{
  vtkIntArray * data = vtkIntArray::New();
  data -> SetNumberOfComponents(1);
  data -> SetName(dataname.c_str());

  for(unsigned int ii=0; ii<cldata.size(); ++ii)
    data -> InsertComponent(ii, 0, cldata[ii]);

  grid_w -> GetCellData() -> AddArray( data );
  data -> Delete();
}

// EOF
