/*
 * file use_metis.h
 * brief convert mesh to graph
 * author lxy
 * date 20170830
*/

#ifndef _USE_METIS_H
#define _USE_METIS_H

#define NODEMAXIUMSIZE 1000000

#include <iterator>
#include <cstring>
#include <string>
#include <limits>
#include <math.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iostream>
#include "metis.h"
#include "surface_function.h"



struct Vector_3d
{
    double x,y,z;
    Vector_3d(double x,double y,double z):x(x),y(y),z(z){}
    Vector_3d() {}
    Vector_3d operator + (const Vector_3d &a) const{return Vector_3d(x+a.x,y+a.y,z+a.z);}
    Vector_3d operator - (const Vector_3d & a) const{ return Vector_3d(x - a.x, y - a.y, z - a.z); }
    Vector_3d operator * (const double & a) const{ return Vector_3d(x * a, y * a, z * a); }
    Vector_3d operator / (const double & a) const{ return Vector_3d(x / a, y / a, z / a); }
    double sqrlen() const{ return x*x + y*y + z*z; }
    double length() const{ return sqrt(sqrlen()); }

};
void neighbors(std::string filename ,std::vector<int> & center, int n_number );
void circle_number(std::string filename, std::vector<int> & center, int n_number ,double zoneRadius);
void layer_number(std::string filename ,std::vector<int> & layer, int n_number );

void mesh2dual (std::string filename, std::vector<int> &ary_xadj, std::vector<int> &ary_adjncy);
void geof_top_bool(std::string filename, std::vector<int> &layer);

void geof2CSR(std::string filename, std::vector<int> &xadj, std::vector<int> &adjncy);
void geof_center_bool(std::string filename, std::vector<int> &center);
void FindCenterNode(std::string fileline, double dist_lxy[]);
void ReadElementDistance(std::string fileline, int &center_number, double &avg_dist, double dist_lxy[]);
void InputLineElement(std::string fileline,std::vector<int> &xadj, std::vector<int> &adjncy);
void FindTopNode(std::string fileline, bool top[]);
void ReadElementNode(std::string fileline, std::vector<int> &layer, bool top[]);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




void InputLineElement(std::string fileline, std::vector<int> &xadj, std::vector<int> &adjncy)
{
    std::istringstream temp;
    temp.str(fileline);
    int n;
    int sum=xadj.back();
    std::string t;
    temp>>n>>t;
    while(temp>>n) {adjncy.push_back(n);sum++;}
    xadj.push_back(sum);

}

void geof2CSR(std::string filename, std::vector<int> &xadj, std::vector<int> &adjncy)
{
    char * cstr = new char [filename.length()+1];
    std::strcpy (cstr, filename.c_str());

    std::ifstream FIn;
    FIn.open(cstr);
    if (! FIn.is_open())
       { std::cout << "Error opening file"; exit (1); }

    std::string fileline;
    xadj.push_back(0);

    bool beginflag=true;
    bool endflag=true;
    while (beginflag)
    {
        getline(FIn,fileline);
        if(fileline.find("**element")!= std::string::npos) {beginflag=false;getline(FIn,fileline);break;}
    }
    while (endflag)
    {
        getline(FIn,fileline);
        if ( (fileline.find("**group")!=std::string::npos) ) {endflag=false;break;}
        InputLineElement(fileline,xadj,adjncy);
    }
    FIn.close();
}


void mesh2dual(std::string filename, std::vector<int> &ary_xadj, std::vector<int> &ary_adjncy)
{


    std::vector<int> axadj,aadjncy;

    geof2CSR(filename,axadj,aadjncy);

    int NN=axadj[axadj.size()-1];
    idx_t ne=axadj.size()-1;

    std::vector<int>::iterator biggest = std::max_element(aadjncy.begin(), aadjncy.end());
    //idx_t nn=*biggest; //以前都用的很正常，现在放到服务器上又说不对，好不容易做的找最大值，原来它要的就是单纯的点被数了多少次，也不知道以前是怎么运行的
    idx_t nn=NN;
    idx_t ncommon=1;
    idx_t numflag=0;

    idx_t eptr[ne+1];
    idx_t eind[NN];

    for (int i=0;i<ne+1;i++) eptr[i]=axadj[i];
    for (int i=0;i<NN;i++) eind[i]=aadjncy[i];

    //idx_t *xadj=std::nullptr;
    //idx_t *adjncy=std::nullptr;
    idx_t *xadj;
    idx_t *adjncy;


    METIS_MeshToDual(&ne,&nn,eptr,eind,&ncommon,&numflag,&xadj,&adjncy);

    int i,j;
    i=0;
    while ( i<ne+1)
    {ary_xadj.push_back(xadj[i]);i++;}// ary_xadj's size = number of elements +1


    j=0;
    while(j<xadj[ne])
    {ary_adjncy.push_back(adjncy[j]);j++;}// ary_adjncy's size = maxium value in xadj

    METIS_Free(xadj);
    METIS_Free(adjncy);
}






void geof_top_bool(std::string filename, std::vector<int> &layer)
{
    char * cstr = new char [filename.length()+1];
    std::strcpy (cstr, filename.c_str());

    std::ifstream FIn;
    FIn.open(cstr);
    if (! FIn.is_open())
       { std::cout << "Error opening file"; exit (1); }

    std::string fileline;
    bool top[1000000];

    bool beginflag=true;
    bool endflag=true;
    bool flag3=true;
    while (beginflag)
    {
        getline(FIn,fileline);
        if(fileline.find("**node")!= std::string::npos) {beginflag=false;getline(FIn,fileline);break;}
    }
    while (endflag)
    {
        getline(FIn,fileline);
        if ( (fileline.find("**element")!=std::string::npos) ) {endflag=false;break;}
        FindTopNode(fileline,top);
    }
    getline(FIn,fileline);
    while (flag3)
    {
        getline(FIn,fileline);
        if ( (fileline.find("**group")!=std::string::npos) ) {flag3=false;break;}
        ReadElementNode(fileline,layer,top);
    }



    FIn.close();
}

void FindTopNode(std::string fileline, bool top[])
{
    std::istringstream temp;
    temp.str(fileline);
    int n;
    double x,y,z;
    std::string t;
    temp>>n>>x>>y>>z;
    top[n-1]=true;
    if ( fabs( y-surface_function(x) )< err_lxy ) {top[n-1]=false;}

}

void ReadElementNode(std::string fileline, std::vector<int> &layer, bool top[])
{
    std::istringstream temp;
    temp.str(fileline);
    int n,noden;
    std::string t;
    temp>>n>>t;
    while(temp>>noden)
    {
        if(top[noden-1]==false) {layer.push_back(n-1);break;}
    }

}




void geof_center_bool(std::string filename, std::vector<int> &center ,std::vector<double> &distance2Center)
{
    char * cstr = new char [filename.length()+1];
    std::strcpy (cstr, filename.c_str());

    std::ifstream FIn;
    FIn.open(cstr);
    if (! FIn.is_open())
       { std::cout << "Error opening file"; exit (1); }

    std::string fileline;
    int center_number;
    double avg_dist;

    int RealCenter;
    double min_dist=(std::numeric_limits<double>::max)();
    bool beginflag=true;
    bool endflag=true;
    bool flag3=true;
    double dist_lxy[NODEMAXIUMSIZE];

    while (beginflag)
    {
        getline(FIn,fileline);
        if(fileline.find("**node")!= std::string::npos) {beginflag=false;getline(FIn,fileline);break;}
    }
    while (endflag)
    {
        getline(FIn,fileline);
        if ( (fileline.find("**element")!=std::string::npos) ) {endflag=false;break;}
        FindCenterNode(fileline,dist_lxy);
    }
    getline(FIn,fileline);
    while (flag3)
    {
        getline(FIn,fileline);
        if ( (fileline.find("**group")!=std::string::npos) ) {flag3=false;break;}
        ReadElementDistance(fileline,center_number,avg_dist,dist_lxy);
        distance2Center.push_back(avg_dist);
        if(avg_dist<min_dist) {min_dist=avg_dist; RealCenter=center_number; }
    }

    center.push_back(RealCenter);


    FIn.close();
}

void FindCenterNode(std::string fileline, double dist_lxy[])
{
    std::istringstream temp;
    temp.str(fileline);
    int n;
    double x,y,z;
    std::string t;
    temp>>n>>x>>y>>z;

    dist_lxy[n-1]=distance_fun(x,y,z);

}
void ReadElementDistance(std::string fileline, int &center_number, double &avg_dist, double dist_lxy[])
{
    std::istringstream temp;
    temp.str(fileline);
    int n,noden,taille;
    std::string t;
    temp>>n>>t;
    taille=0;
    avg_dist=0.0;
    center_number=n-1;
    while(temp>>noden)
    {
        taille++;
        avg_dist+=dist_lxy[noden-1];
    }
    avg_dist=avg_dist/taille;
}




void layer_number(std::string filename ,std::vector<int> & layer, int n_number )
{

    std::vector<int> xadj; std::vector<int> adjncy;
    mesh2dual(filename,xadj,adjncy);
    geof_top_bool(filename,layer);


    int begin_layer,end_layer;
    begin_layer=0;
    end_layer=layer.size();

    while (n_number>1)
    {
        for (int i=begin_layer;i<end_layer;i++)
        {
            int t1,t2;
            t1=xadj[ layer[i] ];t2=xadj[ layer[i] +1]-1;
            for( int i=t1;i<t2;i++)
            { layer.push_back(adjncy[i]); }
        }
        sort(layer.begin(),layer.end());
        layer.erase(unique(layer.begin(),layer.end()),layer.end());

        begin_layer=end_layer;
        end_layer=layer.size();
        n_number--;
    }


}


void circle_number(std::string filename ,std::vector<int> & center, int n_number, double zoneRadius )
{

    std::vector<double> distance2Center;
    std::vector<int> xadj; std::vector<int> adjncy;
    mesh2dual(filename,xadj,adjncy);
    geof_center_bool(filename,center,distance2Center);

    if (distance2Center.size()!= (xadj.size()-1) ) exit(1);

    int begin_circle,end_circle;


    while (n_number>1)
    {
        begin_circle=0;
        end_circle=center.size();
        for (int i=begin_circle;i<end_circle;i++)
        {
            int t1,t2;
            t1=xadj[center[i]];t2=xadj[center[i]+1]-1;
            for( int i=t1;i<t2;i++)
            {
                if ( (distance2Center[ adjncy[i] ])<zoneRadius)   center.push_back(adjncy[i]);

            }
        }
        sort(center.begin(),center.end());
        center.erase(unique(center.begin(),center.end()),center.end());

        begin_circle=end_circle;
        end_circle=center.size();
        n_number--;
    }


}

void neighbors(std::string filename ,std::vector<int> & center, int n_number )
{

    std::vector<double> distance2Center;
    std::vector<int> xadj; std::vector<int> adjncy;
    mesh2dual(filename,xadj,adjncy);

    int begin_circle,end_circle;

    while (n_number>1)
    {
        begin_circle=0;
        end_circle=center.size();
        for (int i=begin_circle;i<end_circle;i++)
        {
            int t1,t2;
            t1=xadj[center[i]];t2=xadj[center[i]+1]-1;
            for( int i=t1;i<t2;i++)
            {
                center.push_back(adjncy[i]);
            }
        }
        sort(center.begin(),center.end());
        center.erase(unique(center.begin(),center.end()),center.end());
        n_number--;
    }


}


#endif // _USE_METIS_H
