#include "headers/fatAnalysis.h"

enum modelType
	{
		isotropic,
		cubic,
		j2,
		crystal
	};
int setType(const int n, modelType &m)
{
	switch(n)
	{
		
	case 1:
		m=isotropic;
		break;
	case 2:
		m=cubic;
		break;
	case 3:
		m=j2;
		break;
	case 4:
		m=crystal;
		break;
	default:
		std::cout<<"Not Correct Type"<<std::endl;
		return -1;
	}
	return 0;
}

int main(int argc,char *argv[])
{
	
	modelType modeltype;
	
	if(argc==1)
	{
		std::cout<<"No Input"<<std::endl;
		return -1;
	}
	
	if(argc==2)
	{
		modeltype=isotropic;
	}
    
    if(argc==3) 
    {
		int temp;
        std::stringstream ss;
        ss<<argv[2];
		ss>>temp;
		ss.clear();
		setType(temp,modeltype);
    }
    
    if(argc==4) 
    {
		int temp;
        std::stringstream ss;
        ss<<argv[2];
		ss>>temp;
		ss.clear();
		setType(temp,modeltype);
        
        ss<<argv[3];
        ss>>MAXIUM_CYCLE;
        ss.clear();
    }
    
    MPI_Init(&argc,&argv);
    
    if(modeltype==isotropic or modeltype==j2)
    {
		HomoModel test;
    
	    double t;
	    double s;
	    
	    t=232.5;
	    s=147.5;
	    
	    double aDV,aPap,aM;
	    aDV=(s-t/2.0)/t*3.0;
	    aPap=3*s/t-std::sqrt(3);
	    aM=2.0*s/t-1.0;
	    
	    test.Init_Elements(argv[1]);
	
	    test.Matake_Analysis(aM,s);
	    test.Papadopoulos_Analysis(aPap,s);    
	    test.DangVan_Analysis(aDV,s); 
	    
	    test.DangVan_Analysis(aDV,s,6);
	    test.Matake_Analysis(aM,s,6);
	    test.Papadopoulos_Analysis(aPap,s,6);
	    
	    
	    test.Fin_Elements();
	}
	
	if(modeltype==cubic)
	{
		PolyModel test;
    
	    double t;
	    double s;
	    
	    t=232.5;
	    s=147.5;
	    
	    double aDV,aPap,aM;
	    aDV=(s-t/2.0)/t*3.0;
	    aPap=3*s/t-std::sqrt(3);
	    aM=2.0*s/t-1.0;
	    
	    test.Init_Grains(argv[1]);
	    test.cubDangVan_Analysis(aDV,s); 
	    test.cubMatake_Analysis(aM,s);
	    test.cubPapadopoulos_Analysis(aPap,s);
	    test.Fin_Grains();
	}
	
	if(modeltype==crystal)
	{
		PolyModel test;
    
	    double t;
	    double s;
	    
	    t=232.5;
	    s=147.5;
	    
	    double aDV,aPap,aM;
	    aDV=(s-t/2.0)/t*3.0;
	    aPap=3*s/t-std::sqrt(3);
	    aM=2.0*s/t-1.0;
	    
	    test.Init_Grains(argv[1]);
	    test.grainDangVan_Analysis(aDV,s); 
	    test.grainMatake_Analysis(aM,s);
	    test.grainPapadopoulos_Analysis(aPap,s);
	    test.Fin_Grains();
		
	}
	
	
	
	
	
	
	
	
	
	
    MPI_Finalize();
    return EXIT_SUCCESS;
}
