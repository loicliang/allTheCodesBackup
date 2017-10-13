#include "use_metis.h"
#include "appendix1.h"
class NumericalResults
{
protected:
    char *filename = new char[20];
    
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    
    int ierr=0;
    int nnodes,nelts,nframes;
    int step,cyc,inc,it;
    double time, time_step, time_cycle;
    int NGrains=0;
    int nesets;
    char *name = new char[20];
    int ninc=0;
public:
    
    int Fin_Analysis()
    {
        VTK_fin();
        dodoF2C_Close();

        return EXIT_SUCCESS;
    }
    
    
    int Init_Analysis(char *file_NAME)
    {
        
      
        std::strcpy(filename,file_NAME);
        dodoF2C_Open(filename,strlen(filename),&ierr);
std::cout<<"***********************"<<std::endl;
std::cout<<"Fatigue Analysis Begins"<<std::endl;
std::cout<<"***********************"<<std::endl;
    
        dodoF2C_GetNodesNumber(&nnodes);
        dodoF2C_GetElementsNumber(&nelts);
        dodoF2C_GetFramesNumber(&nframes);
        dodoF2C_GetEsetsNumber(&nesets);
    
std::cout<<"Nodes:"<<nnodes<<"Elements:"<<nelts<<"Sets:"<<nesets<<"Frames:"<<nframes<<std::endl;  

        ierr=0;
        for (int f=0;f<nframes+1;f++)
        {
            dodoF2C_ReadFrame(f,&ierr);
            if (ierr!=0) {break;}
            dodoF2C_GetIncrement(&inc);
            dodoF2C_GetStep(&step);
            dodoF2C_GetCycle(&cyc);
            dodoF2C_GetTotalTime(&time);
            if (cyc==MAXIUM_CYCLE) {ninc++;}
        }
 
std::cout<<"N_INC"<<ninc<<std::endl; 
    VTK_init();
    return EXIT_SUCCESS;
    }
    
    
    int VTK_init()
    {
        GetMesh(unstructuredGrid,nnodes,nelts);
        return EXIT_SUCCESS;
    }
    
    int VTK_fin()
    {
        std::ostringstream sstream;
        sstream << filename << ".vtu";
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
        
        return EXIT_SUCCESS;
    }    
    
    
};

class HomoModel: private NumericalResults //only elements
{
protected:
   // ELEMENT_lxy element_lxy[nelts];
    ELEMENT_lxy *element_lxy;
public:
    int Fin_Elements()
    {        
        delete [] element_lxy;
        Fin_Analysis();
        return EXIT_SUCCESS;
    }
    
    int Set_ElementSize(int num)
    {
        ELEMENT_lxy *temp= new ELEMENT_lxy[num];
        element_lxy=temp;
        
        return EXIT_SUCCESS;
    }
    
    
    int Init_Elements(char *file_NAME)
    {
        Init_Analysis(file_NAME);
        
        Set_ElementSize(nelts);
        
        
//read volume        
        for (int i=0;i<nelts;i++)  
        {
            element_lxy[i].gr=0;
            double voltemp;
            dodoF2C_GetElementVolume(i+1, &voltemp);
            element_lxy[i].volume=voltemp;
        }
//read stress
        ninc=0;ierr=0;
        for (int f=0;f<nframes+1;f++)
        {
            dodoF2C_ReadFrame(f,&ierr);
            if (ierr!=0) {break;}
            dodoF2C_GetIncrement(&inc);
            dodoF2C_GetStep(&step);
            dodoF2C_GetCycle(&cyc);
            dodoF2C_GetTotalTime(&time);
            if (cyc==MAXIUM_CYCLE) 
            {
                for (int i=0;i<nelts;i++) 
                {
                    if (element_lxy[i].gr!=-1) 
                    { 
                        double tens[6]={0,0,0,0,0,0};
                        tensor2 tens2;
                        dodoF2C_GetElementValueTensor2s(i+1,"sig",3,tens);
                        tens2=tens;
                        element_lxy[i].sig.push_back(tens2);
                        
                    }
                }
            ninc++;
            }
        }
    return EXIT_SUCCESS;
    }
    
    
    int DangVan_Analysis(double alphaDV, double sbeta)
    {
        double dv=0.0;
        for (int i=0;i<nelts;i++)
        {
            DangVan_macro_lxy(ninc,element_lxy[i],alphaDV);
            if (element_lxy[i].dangvan>dv) {dv= element_lxy[i].dangvan;}
        }
        std::cout << "DV " << dv << std::endl;
        
        {
        vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
        vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
        vals->SetName("Dang Van HomoModel"); //set the name of the value
        double scal=0.0;
        for (int i=0;i<nelts;i++)
        {
            scal=(element_lxy[i].dangvan)/sbeta;
            vals->InsertNextValue(scal);          
        }
        unstructuredGrid->GetCellData()->AddArray(vals);
        }
        return EXIT_SUCCESS;
    }
    
    int HuyenMorel_Analysis(double alphaHM, double gamma, double tau0p, int m)
    {
        double huyen=0.0;
        HuyenMorel_macro_lxy (ninc,nelts,element_lxy, alphaHM, gamma, tau0p, m, &huyen);
        std::cout<<"HM "<<huyen<<std::endl;
        {
        vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
        vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
        vals->SetName("Huyen Morel"); //set the name of the value
        double scal=0.0;
        for (int i=0;i<nelts;i++) 
        {
            scal=element_lxy[i].huyen;
            vals->InsertNextValue(scal);
        }
        unstructuredGrid->GetCellData()->AddArray(vals);
        }
        
        return EXIT_SUCCESS;
    }

    int Matake_Analysis(double alphaM, double sbeta)
    {
        double matake=0.0;
        for (int i=0;i<nelts;i++)
        {
            Matake_macro_lxy(ninc,element_lxy[i],alphaM);
            if (element_lxy[i].matake>matake) {matake= element_lxy[i].matake;}
        }
        std::cout << "Matake " << matake << std::endl;
        
        {
        vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
        vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
        vals->SetName("Matake HomoModel"); //set the name of the value
        double scal=0.0;
        for (int i=0;i<nelts;i++)
        {
            scal=(element_lxy[i].matake)/sbeta;
            vals->InsertNextValue(scal);          
        }
        unstructuredGrid->GetCellData()->AddArray(vals);
        }
        return EXIT_SUCCESS;
    }
    
    int Papadopoulos_Analysis(double alphaP,double sbeta)
    {
       double papa=0.0;
       for(int i=0;i<nelts;i++)
       {
           Papadopoulos_macro_lxy(ninc,element_lxy[i],alphaP);
           if (element_lxy[i].papadopoulos > papa) {papa= element_lxy[i].papadopoulos;}
       }
       std::cout << "Papadopoulos " << papa << std::endl;
       {
            vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
            vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
            vals->SetName("Papadoupolos HomoModel"); //set the name of the value
            double scal=0.0;
            for (int i=0;i<nelts;i++)
            {
                scal=(element_lxy[i].papadopoulos)/sbeta;
                vals->InsertNextValue(scal);          
            }
            unstructuredGrid->GetCellData()->AddArray(vals);
           
       }
       return EXIT_SUCCESS; 
    }
    
    
    
    int NonLocalCircleAnalysis()
    {
     
        std::cout<<"***********************************"<<std::endl;
        std::cout<<"   Non Local Analysis Conducting   "<<std::endl;
        std::cout<<"***********************************"<<std::endl;
        std::vector<int> center;
        std::string filename;
        filename="tt.geof";
    
        circle_number(filename, center,DEPTH_SEARCH,ZONE_RADIUS );
        int centersize=center.size();
        ELEMENT_lxy CircleZone[centersize];
        
        for(int i=0;i<centersize;i++)
        {
            CircleZone[i]=element_lxy[center[i]];
        }
        
        for (int certainINC=ninc;certainINC>0;certainINC--)
        {   
            double dangvanC=0;
            double volumeC=0;
            for( int i=0;i<centersize;i++)
            {
                int eleLabel=center[i];
                DangVan_macro_lxy(certainINC,element_lxy[eleLabel],alphaDV);
                dangvanC+=element_lxy[eleLabel].dangvan*element_lxy[eleLabel].volume;
                volumeC+=element_lxy[eleLabel].volume;
            }
            dangvanC=dangvanC/volumeC;
            std::cout <<"at the "<<certainINC<<" INC we have " << "DangVan Circle :" << dangvanC << std::endl;
            if (dangvanC<t) break;    
        }
        
        return EXIT_SUCCESS; 
    }
    
};


class PolyModel: private HomoModel //contains grains
{
protected:

    
    
public:
    
    
    
    
    
};



