/*
 * author: Camille ROBERT, LIANG Xiaoyu 
 * Date: 2017.11.20
 * Graphic output of dodoFEM
 * 
*/

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


#include "mpi.h"










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







/*
extern"C" {
   void dodo2vtk_Init(char*, int, int*, int*, int*);
   void dodo2vtk_GetElementVolume(int num, double *volume);
   void dodo2vtk_GetElementType(int num, int *vtktype);
   void dodo2vtk_GetElementKonnec(int num, int vtkkonec[]);
   void dodo2vtk_GetNodeValueVector(int, char*, int, double vect[3]);
   void dodo2vtk_GetElementValueTensor2s(int num,char* cval,int lcval,double tens[6]);
   void dodo2vtk_GetElementMaterialValueScalar(int num,char* cval,int lcval,char* cmat,int lcmat,double* scal);
   void dodo2vtk_GetElementMaterialValueTensor1(int num,char* cval,int lcval,char* cmat,int lcmat,double tens[3]);
   void dodo2vtk_ReadNextFrame(int* step, int* cyc, int* inc, double* time, double* time_step, double* time_cycle, int* ierr);
   void dodo2vtk_End();
   void dodo2vtk_Rewind();
}
*/

#include "dodoF2C.hpp"






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

    int etype,enodes;
    int konnec[100]; 
    for (int i=0; i<nelts;i++)
    {
     dodoF2C_GetElementTypeVTK(i+1,&etype,&enodes);
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




int AvgGetElementValueTensor2s(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nelts,int nesets, std::string name, std::string newname)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());
    char *nameGr= new char[20];

    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(6); //we will have only 1 value associated with the triangle
    vals->SetName(newname.c_str()); //set the name of the value
    
    double tens[6]={0,0,0,0,0,0};
    double tensGr[nesets][6];
    
    bool bbb[nelts];
        int grNo[nelts];
        double volElt[nelts];
        double volGr[nesets];
        int NGrains=0;
        for (int i=0;i<nelts;i++)  {grNo[i]=99999999;}
        for (int n=0;n<nesets;n++)
        {
            dodoF2C_GetEsetName(n+1,nameGr,20);
            if (strstr(nameGr, "grain") != NULL)
            {
                dodoF2C_GetElementsSet(nameGr,20,bbb);
                for (int e=0;e<nelts;e++) 
                {
                    if (bbb[e]!=0) {grNo[e]=NGrains;}
                }
            NGrains++;//ATTENTION!!! the sequence may not be right for hybrid model (J2+Grains) it is a legacy problem !NONONO, IT'S OK
            }
        }
        
        if (NGrains==0)//ATTENTION!!! no grains found,just quit
        {
            std::cout<<"WRONG!No Grains"<<std::endl;
            return -1;
        }

        
        
        for (int i=0;i<nesets;i++) 
        {
            volGr[i]=0.0;
            for (int k=0;k<6;k++){   tensGr[i][k]=0.0;}
        }
		for (int i=0;i<nelts;i++) 
		{
			double voltemp;
			dodoF2C_GetElementVolume(i+1, &voltemp);
			volElt[i]=voltemp;
			if (grNo[i]!=99999999) {volGr[grNo[i]]+=volElt[i];}
		}
		
		double voltotal=0.0;
		for (int i=0;i<NGrains;i++) {voltotal+=volGr[i];} std::cout<<"Total Volume: "<<voltotal<<std::endl;
				
				
                for (int i=0;i<nelts;i++) 
				{
					dodoF2C_GetElementValueTensor2s(i+1,cname,name.size(),tens);
					if (grNo[i]!=99999999) 
					{ 
						int nt;
						nt=grNo[i];
						
						for (int k=0;k<6;k++)
						{
							tensGr[nt][k]+=tens[k]*volElt[i]/volGr[nt];
						}
						
					}
                
				}
    
    
    
        for (int i=0;i<nelts;i++) 
				{
					tens[0]=0.0;
                    tens[1]=0.0;
                    tens[2]=0.0;
                    tens[3]=0.0;
                    tens[4]=0.0;
                    tens[5]=0.0;
                    
					if (grNo[i]!=99999999) 
					{ 
						int nt;
						nt=grNo[i];
						
						for (int k=0;k<6;k++)
						{
							tens[k]=tensGr[nt][k];
						}
						
					}
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





int main(int argc, char *argv[])
{


 int ierr = MPI_Init ( &argc, &argv );



  if(argc != 2)
  {
   std::cout << "Required arguments: OutputFilename" << std::endl;
   return EXIT_FAILURE;
    MPI_Finalize ( );
  }

   int nnodes,nelts,nframes,nesets;


   dodoF2C_Open(argv[1],strlen(argv[1]),&ierr);

   dodoF2C_GetNodesNumber(&nnodes);
   dodoF2C_GetElementsNumber(&nelts);
   dodoF2C_GetFramesNumber(&nframes);
	dodoF2C_GetEsetsNumber(&nesets);

   std::cout << nnodes << " " << nelts << " " << nframes<<" "<<nesets << std::endl;


   ofstream pvdfile;
   std::string pvdfilename;
   {
    std::ostringstream sstream;
    sstream << argv[1] << ".pvd";
    std::string filepvdname=sstream.str();
    pvdfile.open(filepvdname.c_str());
   }


   pvdfile << "<?xml version=\"1.0\"?>" << std::endl;
   pvdfile << "<VTKFile type=\"Collection\" version=\"0.1\"" << std::endl;
   pvdfile << "         byte_order=\"LittleEndian\"" << std::endl;
   pvdfile << "         compressor=\"vtkZLibDataCompressor\">" << std::endl;
   pvdfile << "  <Collection>" << std::endl;



   int step,cyc,inc,it;
   double time, time_step, time_cycle;
   it=0;
   for (int f=0;f<nframes+1;f++)
   {
    dodoF2C_ReadFrame(f,&ierr);
    if (ierr!=0) {break;}
    //if (ierr!=0) break;
    dodoF2C_GetIncrement(&inc);
    dodoF2C_GetStep(&step);
    dodoF2C_GetCycle(&cyc);
    dodoF2C_GetTotalTime(&time);
    std::cout << f << " " << time << std::endl;
    dodoF2C_GetStepTime(&time_step);
    dodoF2C_GetCycleTime(&time_cycle);

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    GetMesh(unstructuredGrid,nnodes,nelts);

    GetNodeValueVector(unstructuredGrid,nnodes,"U","U");
    GetNodeValueVector(unstructuredGrid,nnodes,"F","F");

    GetElementValueTensor2s(unstructuredGrid,nelts,"sig","sig");
    GetElementValueTensor2s(unstructuredGrid,nelts,"eto","eto");
    GetElementValueTensor2s(unstructuredGrid,nelts,"ein","ein");
    
    
    AvgGetElementValueTensor2s(unstructuredGrid,nelts,nesets,"sig","sigAvg");
    AvgGetElementValueTensor2s(unstructuredGrid,nelts,nesets,"eto","etoAvg");
    AvgGetElementValueTensor2s(unstructuredGrid,nelts,nesets,"ein","einAvg");

/*
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui1", "tau1");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui2", "tau2");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui3", "tau3");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui4", "tau4");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui5", "tau5");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui6", "tau6");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui7", "tau7");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui8", "tau8");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui9", "tau9");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui10", "tau10");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui11", "tau11");
    GetElementMaterialValueScalar(unstructuredGrid,nelts,"Mat4","taui12", "tau12");


    GetElementMaterialValueTensor1(unstructuredGrid,nelts,"Mat4","ncs1", "n1");
*/

/*
    {
    vtkDoubleArray *t = vtkDoubleArray::New();
    t->SetName("TIME");
    t->SetNumberOfTuples(1);
    t->SetTuple1(0, time);
    unstructuredGrid->GetFieldData()->AddArray(t);
    }

    {
    vtkIntArray *c = vtkIntArray::New();
    c->SetName("CYCLE");
    c->SetNumberOfTuples(1);
    c->SetTuple1(0, cyc);
    unstructuredGrid->GetFieldData()->AddArray(c);
    }

    {
    vtkIntArray *c = vtkIntArray::New();
    c->SetName("STEP");
    c->SetNumberOfTuples(1);
    c->SetTuple1(0, step);
    unstructuredGrid->GetFieldData()->AddArray(c);
    }
*/


    //it++;
    {
     std::ostringstream sstream;
     sstream << argv[1] << "." << setfill('0') << setw(5) << f << ".vtu";
     std::string vtufilename= sstream.str();
     std::cout << vtufilename << std::endl;
     vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
     writer->SetFileName(vtufilename.c_str());
     #if VTK_MAJOR_VERSION <= 5
      writer->SetInput(unstructuredGrid);
     #else
      writer->SetInputData(unstructuredGrid);
     #endif
     writer->Write(); 

     //pvdfile << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << vtufilename << "\"/>" << std::endl;


     //std::string vtufilename2 = vtufilename.substr(0, vtufilename.find_last_of("\\/"));


  size_t found=vtufilename.find_last_of("/\\");
  //cout << " folder: " << vtufilename.substr(0,found) << endl;
     //cout << vtufilename.substr(0, vtufilename.find_last_of("\\/")) << endl;
     std::string vtufilename2 = vtufilename.substr(found+1);

     pvdfile << "<DataSet timestep=\"" << std::scientific << time << "\" group=\"\" part=\"0\" file=\"" << vtufilename2 << "\"/>" << std::endl;

    }

   }

   pvdfile << "  </Collection>" << std::endl;
   pvdfile << "</VTKFile>" << std::endl;



  dodoF2C_Close();
  MPI_Finalize ( );
  return EXIT_SUCCESS;


}
