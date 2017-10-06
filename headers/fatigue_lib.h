/*
 *file fatigue_lib
 *brief include all the library of the program Fatigue, declare the essential functions
 *author CRobert, lxy
 *date 20170831
 *
*/

#ifndef _FATIGUE_LIB
#define _FATIGUE_LIB


#include "mpi.h"




#include <vtkVersion.h>
#include <vtkSmartPointer.h>

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include "vtkNew.h"

#include <sstream>
#include <iostream>
#include <sstream>
#include <ctype.h>
#include <string.h>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>


#include <ctype.h>
#include <string.h>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <utility>


#include <sys/stat.h>
#include <vector>      // for NodeDataInfoTemp & CellDataInfoTemp
#include <map>         // for TimeStepValuesMap & NumberOfPolygonsMap
#include <string>      // for TimeStepValuesMap & NumberOfPolygonsMap
 
 
 
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator> 

#include <cstring>
#include <string>


#include "dodoF2C.hpp"
#include "DV.hpp"
#include "Papadopoulos.hpp"
#include "HuyenMorel.hpp"



int GetMesh(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nnodes, int nelts);
int GetNodeValueVector(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nnodes, std::string name, std::string newname);
int GetElementValueTensor2s(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nelts, std::string name, std::string newname);
int GetElementMaterialValueScalar(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nelts, std::string mat, std::string name, std::string newname);
int GetElementMaterialValueTensor1(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nelts, std::string mat, std::string name, std::string newname);



int GetMesh(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nnodes, int nelts)
{

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    std::string name="COO";
    char *cname = new char[name.length() + 1];
    std::strcpy(cname, name.c_str());

    //std::cout << name.length() << " " << cname << std::endl;

    for (int i=0; i<nnodes; i++)
    {
     double vect[3]={0,0,0};
     dodoF2C_GetNodeValueVector(i+1,cname,name.length(),vect);
     //std::cout << vx << " " << vy << " " << vz << std::endl;
     points->InsertNextPoint(vect);
    }
    unstructuredGrid->SetPoints(points);
    delete[] cname;

    int etype;
    int konnec[100]; 
    for (int i=0; i<nelts;i++)
    {
     dodoF2C_GetElementTypeVTK(i+1,&etype);
     dodoF2C_GetElementKonnecVTK(i+1,konnec);

     switch (etype)
     {
      case 5: {
       vtkIdType ptIds[] = {konnec[0]-1, konnec[1]-1, konnec[2]-1};
       unstructuredGrid->InsertNextCell( VTK_TRIANGLE, 3, ptIds );
       break;
      }
      case 10: {
       vtkIdType ptIds[] = {konnec[0]-1, konnec[1]-1, konnec[2]-1, konnec[3]-1};
       unstructuredGrid->InsertNextCell( VTK_TETRA, 4, ptIds );
       break;
      }
      case 12: {
       vtkIdType ptIds[] = {konnec[0]-1, konnec[1]-1, konnec[2]-1, konnec[3]-1, konnec[4]-1, konnec[5]-1, konnec[6]-1, konnec[7]-1};
       unstructuredGrid->InsertNextCell( VTK_HEXAHEDRON, 8, ptIds );
       break;
      }
      case 24: {
       vtkIdType ptIds[] = {konnec[0]-1, konnec[1]-1, konnec[2]-1, konnec[9]-1, konnec[3]-1, konnec[4]-1, konnec[5]-1, konnec[6]-1, konnec[7]-1, konnec[8]-1};
       unstructuredGrid->InsertNextCell( VTK_QUADRATIC_TETRA, 10, ptIds );
       break;
      }
      default: {
       std::cerr << "ERROR, ELEMENT TYPE " << etype << "NOT YET IMPLEMENTED";
      }
     }
    }



  return 0;
}


int GetNodeValueVector(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nnodes, std::string name, std::string newname)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());
    {
     vtkSmartPointer<vtkDoubleArray> vectors =vtkSmartPointer<vtkDoubleArray>::New();
     vectors->SetNumberOfComponents(3);
     vectors->SetName(newname.c_str());
     for (int i=0;i<nnodes;i++)
     {
      double vect[3]={0,0,0};
      dodoF2C_GetNodeValueVector(i+1,cname,name.size(),vect);
      //std::cout << i+1 << " " << vect[0] << " " << vect[1] << " " << vect[2] << std::endl;
      vectors->InsertNextTupleValue(vect);
     }
     unstructuredGrid->GetPointData()->AddArray(vectors);
    }
    delete[] cname;

  return 0;
}


int GetElementValueTensor2s(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nelts, std::string name, std::string newname)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());

    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(6); //we will have only 1 value associated with the triangle
    vals->SetName(newname.c_str()); //set the name of the value
    double tens[6]={0,0,0,0,0,0};
    for (int i=0;i<nelts;i++) {
     dodoF2C_GetElementValueTensor2s(i+1,cname,name.size(),tens);
     vals->InsertNextTupleValue(tens);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);

    delete[] cname;
    return 0;

}


int GetElementMaterialValueScalar(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nelts, std::string mat, std::string name, std::string newname)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());
    char *cmat = new char[mat.length() + 1]; // or
    std::strcpy(cmat, mat.c_str());



    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
    vals->SetName(newname.c_str()); //set the name of the value
    double scal=0.0;
    for (int i=0;i<nelts;i++) {
     dodoF2C_GetElementMaterialValueScalar(i+1,cname,name.size(),cmat,mat.size(),&scal);
     //std::cout << scal << std::endl << std::endl;
     vals->InsertNextValue(scal);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);

    delete[] cname;
    delete[] cmat;
    return 0;

}


int GetElementMaterialValueTensor1(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nelts, std::string mat, std::string name, std::string newname)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());
    char *cmat = new char[mat.length() + 1]; // or
    std::strcpy(cmat, mat.c_str());



    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(3); //we will have only 1 value associated with the triangle
    vals->SetName(newname.c_str()); //set the name of the value
    double tens[3]={0,0,0};
    for (int i=0;i<nelts;i++) {
     dodoF2C_GetElementMaterialValueTensor1(i+1,cname,name.size(),cmat,mat.size(),tens);
     //std::cout << scal << std::endl << std::endl;
     vals->InsertNextTupleValue(tens);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);

    delete[] cname;
    delete[] cmat;
    return 0;

}












#endif
