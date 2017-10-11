#include "headers/fatAnalysis.h"

int main(int argc,char *argv[])
{
    MPI_Init(&argc,&argv);
    
    
    HomoModel test;
    
    double t;
    double s;
    
    t=232.5;
    s=147.5;
    
    double aDV,aPap,aM;
    aDV=(s-t/2.0)/t*3.0;
    aPap=3*t/s-std::sqrt(3);
    aM=2.0*t/s-1.0;
    
    test.Init_Elements(argv[1]);
    test.DangVan_Analysis(aDV,s); 
    test.Matake_Analysis(aM,s);
    test.Papadopoulos_Analysis(aPap,s);
    test.Fin_Elements();

    MPI_Finalize();
    return EXIT_SUCCESS;
}
