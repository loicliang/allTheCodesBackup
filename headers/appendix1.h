#ifndef H_app
#define H_app

#include "cal.h"
#include "fatigue_lib.h"

struct NODE_lxy
{
    int n;
    double x,y,z;
};

struct ELEMENT_lxy
{
 int n;
 int n1,n2,n3;
 int gr;
 
 double volume;
 //double norm[3];
 //double dir[3];

 //std::vector <tensor1> norm;
 //std::vector <tensor1> dir;
 
 //double sig[6];
 std::vector <tensor2> sig;
 bool type;//true-> polycrystal false->J2

 
 double huyen;
 double papadopoulos;
 double dangvan;
 NODE_lxy mid_point;
};

struct GRAIN_lxy
{
  int n;
  double volume;
  double frac;
  std::vector <int>  subset;
   //double norm[3];
 //double dir[3];
 std::vector <tensor1> norm;
 std::vector <tensor1> dir;
 
 //double sig[6];
 std::vector <tensor2> sig;
  double huyen;
  double papadopoulos;
  double dangvan;
  NODE_lxy mid_point;  
    
};





void HuyenMorel_macro_lxy (int ninc, int nelts, ELEMENT_lxy element_lxy[], double alpha, double gamma, double tau0p, int m, double *huyen_marco);

void HuyenMorel_polycrystal_lxy (int ninc, int NGrains, int nsys, GRAIN_lxy grain_lxy[], double alpha, double gamma, double tau0p, int m, double *huyen_poly);

void Papadopoulos_polycrystal_lxy(int ninc, int NGrains, int nsys, GRAIN_lxy grain_lxy[], double alpha, double *papadppoulos);

void DangVan_crystal_lxy(int ninc, int nsys, GRAIN_lxy &grain_lxy, double alpha);

void DangVan_macro_lxy(int ninc,ELEMENT_lxy &element_lxy, double alpha);

/* 函数主体 */
void HuyenMorel_macro_lxy (int ninc, int nelts, ELEMENT_lxy element_lxy[], double alpha, double gamma, double tau0p, int m, double *huyen_marco)
{
    int i,j,k;
    double sig[nelts][ninc][6];
    double vol[nelts];
    for (i=0;i<nelts;i++)
    {
        vol[i]=element_lxy[i].volume;
        for (j=0;j<ninc;j++)
        {
            sig[i][j][0]=element_lxy[i].sig[j].a11;
            sig[i][j][1]=element_lxy[i].sig[j].a22;
            sig[i][j][2]=element_lxy[i].sig[j].a33;
            sig[i][j][3]=element_lxy[i].sig[j].a21;
            sig[i][j][4]=element_lxy[i].sig[j].a23;
            sig[i][j][5]=element_lxy[i].sig[j].a13;
        }
    }
    
    double huyenM[nelts];
    
    F2C_HuyenMorel_MACRO(ninc,nelts,sig,vol,alpha,gamma,tau0p,m,huyen_marco,huyenM);
    
    
    for (i=0;i<nelts;i++)
    {
        element_lxy[i].huyen=huyenM[i];
    }
    
    
}


void HuyenMorel_polycrystal_lxy (int ninc, int NGrains, int nsys, GRAIN_lxy grain_lxy[], double alpha, double gamma, double tau0p, int m, double *huyen_poly)
{
    int i,j,k;
    double sigG[NGrains][ninc][6];
    double norm[NGrains][nsys][3];
    double dir[NGrains][nsys][3];
    for (i=0;i<NGrains;i++)
    {
        for (j=0;j<ninc;j++)
        {
            sigG[i][j][0]=grain_lxy[i].sig[j].a11;
            sigG[i][j][1]=grain_lxy[i].sig[j].a22;
            sigG[i][j][2]=grain_lxy[i].sig[j].a33;
            sigG[i][j][3]=grain_lxy[i].sig[j].a21;
            sigG[i][j][4]=grain_lxy[i].sig[j].a23;
            sigG[i][j][5]=grain_lxy[i].sig[j].a13;
        }
    }
    
    for (i=0;i<NGrains;i++)
    {
        for (j=0;j<nsys;j++)
        {
            norm[i][j][0]=grain_lxy[i].norm[j].a1;
            norm[i][j][1]=grain_lxy[i].norm[j].a2;
            norm[i][j][2]=grain_lxy[i].norm[j].a3;
         
            dir[i][j][0]=grain_lxy[i].dir[j].a1;
            dir[i][j][1]=grain_lxy[i].dir[j].a2;
            dir[i][j][2]=grain_lxy[i].dir[j].a3;
            
        }
    }
    
    double huyenG[NGrains];
 
    
    F2C_HuyenMorel_POLYCRYSTAL(ninc,NGrains,nsys, sigG, norm, dir, alpha, gamma, tau0p, m, huyen_poly, huyenG);
    
    for (i=0;i<NGrains;i++)
    {
        grain_lxy[i].huyen=huyenG[i];
    }
    
    
}


void Papadopoulos_polycrystal_lxy(int ninc, int NGrains, int nsys, GRAIN_lxy *grain_lxy, double alpha, double *papadopoulos)
{
    int i,j,k;
    double sigG[NGrains][ninc][6];
    double norm[NGrains][nsys][3];
    double dir[NGrains][nsys][3];
    double fracG[NGrains];
    for (i=0;i<NGrains;i++)
    {
        fracG[i]=grain_lxy[i].frac;
        for (j=0;j<ninc;j++)
        {
            sigG[i][j][0]=grain_lxy[i].sig[j].a11;
            sigG[i][j][1]=grain_lxy[i].sig[j].a22;
            sigG[i][j][2]=grain_lxy[i].sig[j].a33;
            sigG[i][j][3]=grain_lxy[i].sig[j].a21;
            sigG[i][j][4]=grain_lxy[i].sig[j].a23;
            sigG[i][j][5]=grain_lxy[i].sig[j].a13;
        }
    }
    
    for (i=0;i<NGrains;i++)
    {
        for (j=0;j<nsys;j++)
        {
            norm[i][j][0]=grain_lxy[i].norm[j].a1;
            norm[i][j][1]=grain_lxy[i].norm[j].a2;
            norm[i][j][2]=grain_lxy[i].norm[j].a3;
         
            dir[i][j][0]=grain_lxy[i].dir[j].a1;
            dir[i][j][1]=grain_lxy[i].dir[j].a2;
            dir[i][j][2]=grain_lxy[i].dir[j].a3;
            
        }
    }
    
    double papadopoulosG[NGrains];
    F2C_Papadopoulos_POLYCRYSTAL(ninc,NGrains,nsys,sigG,norm,dir,fracG,alpha,papadopoulos,papadopoulosG);
    for (i=0;i<NGrains;i++)
    {
        grain_lxy[i].papadopoulos=papadopoulosG[i];
    }
    
    
}

void DangVan_crystal_lxy(int ninc, int nsys, GRAIN_lxy &ggrain_lxy, double alpha)
{
    int i,j,k;
    double sigG[ninc][6];
    double norm[nsys][3];
    double dir[nsys][3];
  
    {
        for (j=0;j<ninc;j++)
        {
            sigG[j][0]=ggrain_lxy.sig[j].a11;
            sigG[j][1]=ggrain_lxy.sig[j].a22;
            sigG[j][2]=ggrain_lxy.sig[j].a33;
            sigG[j][3]=ggrain_lxy.sig[j].a21;
            sigG[j][4]=ggrain_lxy.sig[j].a23;
            sigG[j][5]=ggrain_lxy.sig[j].a13;
        }
    }
    
   
    {
        for (j=0;j<nsys;j++)
        {
            norm[j][0]=ggrain_lxy.norm[j].a1;
            norm[j][1]=ggrain_lxy.norm[j].a2;
            norm[j][2]=ggrain_lxy.norm[j].a3;
         
            dir[j][0]=ggrain_lxy.dir[j].a1;
            dir[j][1]=ggrain_lxy.dir[j].a2;
            dir[j][2]=ggrain_lxy.dir[j].a3;
            
        }
    }
    double dangvanpass;
    F2C_DV_CRYSTAL1(ninc,nsys, sigG, norm, dir, alpha, &dangvanpass);
    ggrain_lxy.dangvan=dangvanpass;
}

void DangVan_macro_lxy(int ninc,ELEMENT_lxy &element_lxy, double alpha)
{
    int i,j,k;
    double sig[ninc][6];

    {
        for (j=0;j<ninc;j++)
        {
            sig[j][0]=element_lxy.sig[j].a11;
            sig[j][1]=element_lxy.sig[j].a22;
            sig[j][2]=element_lxy.sig[j].a33;
            sig[j][3]=element_lxy.sig[j].a21;
            sig[j][4]=element_lxy.sig[j].a23;
            sig[j][5]=element_lxy.sig[j].a13;
        }
    }

    double dangvanpass;
    F2C_DV_Macro(ninc,sig, alpha, &dangvanpass);
    element_lxy.dangvan=dangvanpass;
}





#endif

