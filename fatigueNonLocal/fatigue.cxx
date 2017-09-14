//edited 15:32 20170828



#include "use_metis.h"
#include "appendix1.h"

/*
Main Program
*/

int main(int argc, char *argv[])
{

 int ierr = MPI_Init ( &argc, &argv );



  if(argc != 2)
  {
   std::cout << "Required arguments: OutputFilename" << std::endl;
   return EXIT_FAILURE;
    MPI_Finalize ( );
  }

 

/* Read output file */
   dodoF2C_Open(argv[1],strlen(argv[1]),&ierr);
std::cout<<std::endl;
std::cout<<std::endl;
std::cout<<std::endl;
std::cout<<"***********************"<<std::endl;
std::cout<<"Fatigue Analysis Begins"<<std::endl;
std::cout<<"***********************"<<std::endl;
  
    int nnodes,nelts,nframes;
    dodoF2C_GetNodesNumber(&nnodes);
    dodoF2C_GetElementsNumber(&nelts);
    dodoF2C_GetFramesNumber(&nframes);

   int step,cyc,inc,it;
   double time, time_step, time_cycle;

   int NGrains=0;
   int nesets;
   char *name = new char[20];
   dodoF2C_GetEsetsNumber(&nesets);

std::cout<<"Nodes:"<<nnodes<<"Elements:"<<nelts<<"Sets:"<<nesets<<"Frames:"<<nframes<<std::endl;

//determine the increments involved in calculation
   int ninc=0;
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
   ELEMENT_lxy element_lxy[nelts];
std::cout<<std::endl;    
std::cout<<"elements initialized"<<std::endl;

//find out each grain belongs to which element set
   bool bbb[nelts];
  // int gr[nelts];
   for (int i=0;i<nelts;i++)  {element_lxy[i].gr=-1;}
   
   for (int n=0;n<nesets;n++)
   {
    dodoF2C_GetEsetName(n+1,name, 20);
    if (strstr(name, "grain") != NULL)
    {
     dodoF2C_GetElementsSet(name,20,bbb);
     for (int e=0;e<nelts;e++) 
     {
      //std::cout << bbb[e] << std::endl;
      if (bbb[e]!=0) {element_lxy[e].gr=NGrains;}
     }
     NGrains++;
    }
   }
//if it is only J2 (NGrains==0)
bool J2bool=true;
   if (NGrains==0)
   {
    for (int i=0;i<nelts;i++)  {element_lxy[i].gr=0;}
    NGrains++;
    J2bool=false;
   }
   
   GRAIN_lxy grain_lxy[NGrains];
std::cout<<std::endl;    
std::cout<<"grains initialized"<<std::endl;   
   
   //for (int e=0;e<nelts;e++) {std::cout << gr[e] << std::endl;}
 
//determine the volume of each element and grain
   //double vol[nelts];
   //double volG[NGrains];
   for (int i=0;i<NGrains;i++) {grain_lxy[i].volume=0.0;}
   for (int i=0;i<nelts;i++) 
   {
    double voltemp;
    dodoF2C_GetElementVolume(i+1, &voltemp);
    element_lxy[i].volume=voltemp;
    if (element_lxy[i].gr!=-1) {grain_lxy[element_lxy[i].gr].volume+=element_lxy[i].volume;}
   }

   //double fracG[NGrains];
   
    double voltotal=0.0;
    for (int i=0;i<NGrains;i++) {voltotal+=grain_lxy[i].volume;}
    for (int i=0;i<NGrains;i++) {grain_lxy[i].frac=grain_lxy[i].volume/voltotal;}
std::cout<<std::endl;    
std::cout<<"volumes calculated"<<std::endl;    




int nsys=1;

if (J2bool)
{
   
   {
    int e;
    for (int j=0;j<nelts;j++) {
     if (element_lxy[j].gr>0) {e=j; break;}                       //e？？？？？？？？？？？？？？？？？？？？？？？？？？？？
    }
    dodoF2C_GetElementMaterialValueInteger(e+1,"nsys",4,"Mat4",4,&nsys);
   }
std::cout<<std::endl;    
std::cout<<"slip systems calculated"<<std::endl;    
   //std::cout << "nsys " << nsys << std::endl;


   //double sig[nelts][ninc][6];
   //double sigG[NGrains][ninc][6];
   //double norm[NGrains][nsys][3];
   //double dir[NGrains][nsys][3];

   
    int e;
    for (int i=0;i<NGrains;i++) 
    {
     for (int j=0;j<nelts;j++) 
     {
      if (element_lxy[j].gr==i) {e=j; break;}
     }
     //std::cout << e << std::endl;
     double t1tempa[nsys][3]; tensor1 t1tempt;
     dodoF2C_GetElementMaterialValueVTensor1(e+1,"ncs",3,"Mat4",4,nsys,t1tempa);
    
     for (int k=0;k<nsys;k++)
     {t1tempt=t1tempa[k];grain_lxy[i].norm.push_back(t1tempt);}
     
     
     dodoF2C_GetElementMaterialValueVTensor1(e+1,"bcs",3,"Mat4",4,nsys,t1tempa);
     for (int k=0;k<nsys;k++)
     {t1tempt=t1tempa[k];grain_lxy[i].dir.push_back(t1tempt);}
    }
std::cout<<std::endl;    
std::cout<<"norm,dir calculated"<<std::endl;   
    
}



//    for (int i=0;i<NGrains;i++) 
//    {
//     for (int j=0;j<ninc;j++) 
//     {
//      for (int k=0;k<6;k++) {sigG[i][j][k]=0.0;}
//     }
//    }
//    
//    for (int i=0;i<nelts;i++) 
//    {
//     for (int j=0;j<ninc;j++) 
//     {
//      for (int k=0;k<6;k++) {sig[i][j][k]=0.0;}
//     }
//    }


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
     //bool create=true;   
    
     for (int i=0;i<nelts;i++) 
     {
      //std::cout << gr[i] << std::endl;
      if (element_lxy[i].gr!=-1) 
      { 
        double tens[6]={0,0,0,0,0,0};
        tensor2 tens2;
        tens2=tens;  
        grain_lxy[element_lxy[i].gr].sig.push_back(tens2);
        dodoF2C_GetElementValueTensor2s(i+1,"sig",3,tens);
       //std::cout << tens << std::endl;
       //for (int k=0;k<6;k++) {sigG[gr[i]][ninc][k]+=tens[k]*vol[i]/volG[gr[i]];sig[i][ninc][k]=tens[k];}
        tens2=tens;
        element_lxy[i].sig.push_back(tens2);
        grain_lxy[element_lxy[i].gr].sig[ninc]=grain_lxy[element_lxy[i].gr].sig[ninc]+tens2*element_lxy[i].volume/grain_lxy[element_lxy[i].gr].volume; 
        }
     }
     ninc++;
    }
   }
std::cout<<std::endl;    
std::cout<<"sigma calculated"<<std::endl;
   



//    for (int g=0;g<NGrains;g++)
//    {
//    for (int i=0;i<nsys;i++)
//    {
//    //std::cout << norm[g][i][0] << " " << norm[g][i][1] << " " << norm[g][i][2] << std::endl;
//    }
//    }

   
   double s=232.5;
   double t=147.5;

   vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
   GetMesh(unstructuredGrid,nnodes,nelts);


 
std::cout<<NGrains<<"Grains"<<std::endl; 
std::cout<<std::endl;
std::cout<<std::endl;
std::cout<<std::endl;
std::cout<<"***********************************"<<std::endl;
std::cout<<"ALL Elements and Grains Initialized"<<std::endl;
std::cout<<"***********************************"<<std::endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
if (J2bool)
{
   //double dvG[NGrains];
   {
   double dv=0.0;
   double alpha=(t-s/2.0)/(s/3.0);
   for (int i=0;i<NGrains;i++)
   {
    DangVan_crystal_lxy(ninc,nsys, grain_lxy[i], alpha);
       //F2C_DV_CRYSTAL1(ninc,nsys, sigG[i], norm[i], dir[i], alpha, &dvG[i]);
    if (grain_lxy[i].dangvan>dv) {dv= grain_lxy[i].dangvan;}
    //std::cout << "DV G" << dvG[i] << std::endl;
   }
    std::cout << "DV " << dv << std::endl;
   }
   
      {
    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
    vals->SetName("Dang Van"); //set the name of the value
    double scal=0.0;
    for (int i=0;i<nelts;i++) {
     if (element_lxy[i].gr==-1) {scal=0.0;}
     else {scal=grain_lxy[element_lxy[i].gr].dangvan;}
     vals->InsertNextValue(scal);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);
   }
/*
   double dv[nelts];
   for (int i=0;i<nelts;i++)
   {
    if (gr[i]!=-1) {dv[i]=0.0;}
    else {dv[i]=dvG[gr[i]];}
   }
*/



   //double papadopulosG[NGrains];
   {
   double papadopoulos=0.0;
   double alpha=(3.0*s/t)-std::sqrt(3.0);
   Papadopoulos_polycrystal_lxy( ninc, NGrains, nsys, grain_lxy, alpha, &papadopoulos);
   //F2C_Papadopoulos_POLYCRYSTAL(ninc,NGrains,nsys,sigG,norm,dir,fracG,alpha,&papadopoulos,papadopulosG);
   //for (int i=0;i<NGrains;i++) {std::cout << "PAPADOPOULOS G " << papadopulosG[i] << std::endl;}
    std::cout << "PAPADOPOULOS :" << papadopoulos << std::endl;
   }
     {
    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
    vals->SetName("Papadopoulos"); //set the name of the value
    double scal=0.0;
    for (int i=0;i<nelts;i++) {
     if (element_lxy[i].gr==-1) {scal=0.0;}
     else {scal=grain_lxy[element_lxy[i].gr].papadopoulos;}
     vals->InsertNextValue(scal);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);
   } 
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double huyen=0.0;
   //double huyenG[NGrains];
   //double huyenM[nelts];
   {
   
   //double gamma=2.52e-3;
   //double alpha=0.129;
   //double tau0p=475.0;
   //int m=5;

   double alpha=0.111;
   double gamma=2.45e-3;
   double tau0p=242.0/1.0;
   int m=5;

 if (J2bool)
 {
     HuyenMorel_polycrystal_lxy (ninc,NGrains,nsys,grain_lxy,alpha, gamma, tau0p,m, &huyen);
     //F2C_HuyenMorel_POLYCRYSTAL(ninc,NGrains,nsys, sigG, norm, dir, alpha, gamma, tau0p, m, &huyen, huyenG);
     
        {
    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
    vals->SetName("Huyen Morel"); //set the name of the value
    double scal=0.0;
    for (int i=0;i<nelts;i++) {
     if (element_lxy[i].gr==-1) {scal=0.0;}
     else {scal=grain_lxy[element_lxy[i].gr].huyen;}
     vals->InsertNextValue(scal);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);
   }

}
 else 
 {
    HuyenMorel_macro_lxy (ninc,nelts,element_lxy, alpha, gamma, tau0p, m, &huyen);
     
     //F2C_HuyenMorel_MACRO(ninc,nelts,sig,vol,alpha,gamma,tau0p,m,&huyen,huyenM);
    
     {
    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
    vals->SetName("Huyen Morel"); //set the name of the value
    double scal=0.0;
    for (int i=0;i<nelts;i++) {
     scal=element_lxy[i].huyen;
     
     vals->InsertNextValue(scal);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);
   }
   
   
    {
   double dv=0.0;
   double alpha=(t-s/2.0)/(s/3.0);
   for (int i=0;i<nelts;i++)
   {
    DangVan_macro_lxy(ninc,element_lxy[i],alpha);

    if (element_lxy[i].dangvan>dv) {dv= element_lxy[i].dangvan;}
   }
    std::cout << "DV " << dv << std::endl;
   }
   
    {
    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
    vals->SetName("Dang Van J2"); //set the name of the value
    double scal=0.0;
    for (int i=0;i<nelts;i++)
    {
        scal=element_lxy[i].dangvan;
        vals->InsertNextValue(scal);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);
   }
   
   
   

 }
   
   //for (int i=0;i<NGrains;i++) {std::cout << "HUYEN G " << huyenG[i] << std::endl;}
    std::cout << "HUYEN :" << huyen << std::endl;
   }




//

{
std::cout<<std::endl;
std::cout<<std::endl;
std::cout<<std::endl;
std::cout<<"***********************************"<<std::endl;
std::cout<<"    Selected Analysis Conducting"<<std::endl;
std::cout<<"***********************************"<<std::endl;

//int n_number=5;
std::vector<int>  layer;
std::vector<int> center;
std::string filename;
filename="tt.geof";

circle_number(filename, center,DEPTH_SEARCH,ZONE_RADIUS );
//layer_number(filename ,layer, 1);


std::cout<<"Metis_API excuted"<<std::endl;
//std::vector<ELEMENT_lxy> CircleZone;
//std::vector<ELEMENT_lxy> LayerZone;

int centersize=center.size();
//int layersize=layer.size();

ELEMENT_lxy CircleZone[centersize];
//ELEMENT_lxy LayerZone[layersize];



for(int i=0;i<centersize;i++)
{
    CircleZone[i]=element_lxy[center[i]];
}

//for(int i=0;i<layersize;i++) {    LayerZone[i]=element_lxy[layer[i]];}

std::cout<<"Center Size:   "<<centersize<<std::endl;
//std::cout<<"Layer Size     "<<layersize<<std::endl;
std::cout<<"Zone Selected"<<std::endl;

    double huyenC=0.0;
    double huyenL=0.0;
   //double gamma=2.52e-3;
   //double alpha=0.129;
   //double tau0p=475.0;
   //int m=5;

   double alpha=0.111;
   double gamma=2.45e-3;
   double tau0p=242.0/1.0;
   int m=5;

for (int certainINC=ninc;certainINC>0;certainINC--)
{   
HuyenMorel_macro_lxy (certainINC,centersize,CircleZone, alpha, gamma, tau0p, m, &huyenC);
//HuyenMorel_macro_lxy (certainINC,layersize,LayerZone, alpha, gamma, tau0p, m, &huyenL);
std::cout <<"at the "<<certainINC<<" INC we have " << "HUYEN Circle :" << huyenC << std::endl;
if (huyenC <0.5) break;
//std::cout << "HUYEN Layer: " << huyenL << std::endl;
}

double dangvanC=0;
double volumeC=0;
//double alphaDV=(t-s/2.0)/(s/3.0);  already in cal.h

for (int certainINC=ninc;certainINC>0;certainINC--)
{   
    
    for( int i=0;i<centersize;i++)
    {
        int eleLabel=center[i];
        DangVan_macro_lxy(certainINC,element_lxy[eleLabel],alphaDV);
        dangvanC+=element_lxy[eleLabel].dangvan*element_lxy[eleLabel].volume;
        volumeC+=element_lxy[eleLabel].volume;
    }
    dangvanC=dangvanC/volumeC;
std::cout <<"at the "<<certainINC<<" INC we have " << "DangVan Circle :" << dangvanC << std::endl;
if (dangvanC<s) break;    
    
}
    
    

    
}

//


   //return 0;








    {
     std::ostringstream sstream;
     sstream << argv[1] << ".vtu";
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
    }


  dodoF2C_Close();
  MPI_Finalize ( );
  return EXIT_SUCCESS;


}














