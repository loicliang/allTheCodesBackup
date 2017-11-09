#ifndef H_cal
#define H_cal

#include <string>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <vector>
#include <limits>
//#define alphaDV 0.403226
// #define stabilized_cycle 5
// const double SIGMA=152.5;
struct tensor1;
struct tensor2;
double maxdouble=(std::numeric_limits<double>::max)(); 
double mindouble=(std::numeric_limits<double>::min)();


const int MAXN = 1000;
const double eps = 1e-5;
const double PI = atan2(0.0, -1.0);
inline double sqr(double x){ return x * x; }
inline bool zero(double x){ return (x > 0 ? x : -x) < eps; }
inline int sgn(double x){ return (x > eps ? 1 : (x + eps < 0 ? -1 : 0)); }


struct point3{
    double x, y, z;
    point3(double x, double y, double z):x(x), y(y), z(z){}
    point3() {}
    bool operator == (const point3 & a) const{ return sgn(x - a.x) == 0 && sgn(y - a.y) == 0 && sgn(z - a.z) == 0; }
    bool operator != (const point3 & a) const{ return sgn(x - a.x) != 0 || sgn(y - a.y) != 0 || sgn(z - a.z) != 0; }
    bool operator < (const point3 & a) const{ return sgn(x - a.x) < 0 || (sgn(x - a.x) == 0 && sgn(y - a.y) < 0) || (sgn(x - a.x) == 0 && sgn(y - a.y) == 0 && sgn(z - a.z) < 0); }
    point3 operator + (const point3 & a) const{ return point3(x + a.x, y + a.y, z + a.z); }
    point3 operator - (const point3 & a) const{ return point3(x - a.x, y - a.y, z - a.z); }
    point3 operator * (const double & a) const{ return point3(x * a, y * a, z * a); }
    point3 operator / (const double & a) const{ return point3(x / a, y / a, z / a); }
    point3 operator * (const point3 & a) const{ return point3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x); }//xmult
    double operator ^ (const point3 & a) const{ return x * a.x + y * a.y + z * a.z; }//dmult
    double sqrlen() const{ return sqr(x) + sqr(y) + sqr(z); }
    double length() const{ return sqrt(sqrlen()); }
    point3 trunc(double a) const{ return (*this) * (a / length()); }
    point3 rotate(const point3 & a, const point3 & b, const point3 & c) const{ return point3(a ^ (*this), b ^ (*this), c ^ (*this)); }//abc正交且模为1
};
inline point3 pvec(const point3 & a, const point3 & b, const point3 & c){
    return (a - b) * (b - c);
}//平面法向量
inline bool dotsInline(const point3 & a, const point3 & b, const point3 & c){
    return zero(((a - b) * (b - c)).length());
}//判三点共线
inline bool dotOnlineIn(const point3 & p, const point3 & l1, const point3 & l2){
    return zero(((p - l1) * (p - l2)).length()) && (l1.x - p.x) * (l2.x - p.x) < eps && (l1.y - p.y) * (l2.y - p.y) < eps && (l1.z - p.z) * (l2.z - p.z) < eps;
}//判点是否在线段上,包括端点和共线
inline int decideSide(const point3 & p1, const point3 & p2, const point3 & l1, const point3 & l2){
    return sgn(((l1 - l2) * (p1 - l2)) ^ ((l1 - l2) * (p2 - l2)));
}//点p1和p2,直线l1-l2,-1表示在异侧,0表示在线上,1表示同侧,须保证所有点共面
inline int decideSide(const point3 & p1, const point3 & p2, const point3 & a, const point3 & b, const point3 & c){
    return sgn((pvec(a, b, c) ^ (p1 - a)) * (pvec(a, b, c) ^ (p2 - a)));
}//点p1和p2,平面abc,-1表示在异侧,0表示在面上,1表示同侧
inline bool parallel(const point3 & u1, const point3 & u2, const point3 & v1, const point3 & v2){
    return zero(((u1 - u2) * (v1 - v2)).length());
}//判两直线平行
inline bool parallel(const point3 & a, const point3 & b, const point3 & c, const point3 & d, const point3 & e, const point3 & f){
    return zero((pvec(a, b, c) * pvec(d, e, f)).length());
}//判两平面平行
inline bool parallel(const point3 & l1, const point3 & l2, const point3 & a, const point3 & b, const point3 & c){
    return zero((l1 - l2) ^ pvec(a, b, c));
}//判直线与平面平行
point3 intersection(const point3 & l1, const point3 & l2, const point3 & a, const point3 & b, const point3 & c){
    point3 temp = pvec(a, b, c); return l1 + (l2 - l1) * ((temp ^ (a - l1)) / (temp ^ (l2 - l1)));
}//计算直线与平面交点,须预判是否平行,并保证三点不共线
void intersection(const point3 & a, const point3 & b, const point3 & c, const point3 & d, const point3 & e, const point3 & f, point3 & p1, point3 & p2){
    p1 = parallel(d, e, a, b, c) ? intersection(e, f, a, b, c) : intersection(d, e, a, b, c);
    p2 = parallel(f, d, a, b, c) ? intersection(e, f, a, b, c) : intersection(f, d, a, b, c);
}//计算两平面交线,注意事先判断是否平行,并保证三点不共线,p-q为交线
inline point3 ptoline(const point3 & p, const point3 & l1, const point3 & l2){
    point3 temp = l2 - l1; return l1 + temp * ((p - l1) ^ temp) / (temp ^ temp);
}//点到直线最近点
inline point3 ptoplane(const point3 & p, const point3 & a, const point3 & b, const point3 & c){
    return intersection(p, p + pvec(a, b, c), a, b, c);
}//点到平面最近点
inline double dislinetoline(const point3 & u1, const point3 & u2, const point3 & v1, const point3 & v2){
    point3 temp = (u1 - u2) * (v1 - v2); return fabs((u1 - v1) ^ temp) / temp.length();
}//直线到直线距离
void linetoline(const point3 & u1, const point3 & u2, const point3 & v1, const point3 & v2, point3 & p1, point3 & p2){
    point3 ab = u2 - u1, cd = v2 - v1, ac = v1 - u1;
    p2 = v1 + cd * (((ab ^ cd) * (ac ^ ab) - (ab ^ ab) * (ac ^ cd)) / ((ab ^ ab) * (cd ^ cd) - sqr(ab ^ cd)));
    p1 = ptoline(p2, u1, u2);
}//直线到直线的最近点对,p1在u上,p2在v上,须保证直线不平行



// double sigmah(double t)
// {
//     double temp;
//     t=t-floor(t);
//     if(t>=0.0 && t<=0.25) temp=4.0*t*SIGMA;
//         else if(t>0.25 && t<=0.75) temp=-4.0*SIGMA*(t-0.5);
//         else if(t>0.75 && t<1.0)  temp=4.0*SIGMA*(t-1.0);
//         else temp=0.0;
//     return temp;
// };
// struct array12
// {
//     double a[12];
//     void operator = (const double t[12])
//     {
//         for (int i=0;i<12;i++)
//         {
//             this->a[i]=t[i];
//         }
//     }
// };
// 
// 
// struct Elt
// {
//   int gn=0;
//   int num=0;
//   bool f=false;
//   double volume=0.0;
// };
// struct Grain
// {
//     int num=0;
//     int elts=0;
//     bool f=false;
//     double volume=0.0;
//     std::vector<int> n_elt;
// };



 struct tensor1
{
    double a1,a2,a3;

    double operator * (const tensor1 &A)
    {   double temp;
        temp=this->a1*A.a1+this->a2*A.a2+this->a3*A.a3;
        return temp;
    }
    tensor1 operator * (const double &A)
    {   tensor1 temp;
        temp.a1=this->a1*A;
        temp.a2=this->a2*A;
        temp.a3=this->a3*A;
        return temp;
    }
    tensor1 operator + (const tensor1 &A)
    {   tensor1 temp;
        temp.a1=this->a1+A.a1;
        temp.a2=this->a2+A.a2;
        temp.a3=this->a3+A.a3;
        return temp;
    }
    tensor1 operator - (const tensor1 &A)
    {   tensor1 temp;
        temp.a1=this->a1-A.a1;
        temp.a2=this->a2-A.a2;
        temp.a3=this->a3-A.a3;
        return temp;
    }


    void operator = (const double A[3])
    {

        this->a1=A[0];
        this->a2=A[1];
        this->a3=A[2];
    }

    double mode (void)
    {   double temp;
        temp=sqrt(this->a1*this->a1+this->a2*this->a2+this->a3*this->a3);
        return temp;
    }

};

struct tensor2
{
    double a11,a12,a13,a21,a22,a23,a31,a32,a33;

    tensor1 operator * (const tensor1 &A)
    {   tensor1 temp;
        temp.a1=this->a11*A.a1+this->a12*A.a2+this->a13*A.a3;
        temp.a2=this->a21*A.a1+this->a22*A.a2+this->a23*A.a3;
        temp.a3=this->a31*A.a1+this->a32*A.a2+this->a33*A.a3;
        return temp;
    }
    tensor2 operator + (const tensor2 &A)
    {
        tensor2 temp;
        temp.a11=this->a11+A.a11;
        temp.a12=this->a12+A.a12;
        temp.a13=this->a13+A.a13;
        temp.a21=this->a21+A.a21;
        temp.a22=this->a22+A.a22;
        temp.a23=this->a23+A.a23;
        temp.a31=this->a31+A.a31;
        temp.a32=this->a32+A.a32;
        temp.a33=this->a33+A.a33;
        return temp;
    }
    tensor2 operator / (const double &A)
    {
        tensor2 temp;
        temp.a11=this->a11/A;
        temp.a12=this->a12/A;
        temp.a13=this->a13/A;
        temp.a21=this->a21/A;
        temp.a22=this->a22/A;
        temp.a23=this->a23/A;
        temp.a31=this->a31/A;
        temp.a32=this->a32/A;
        temp.a33=this->a33/A;
        return temp;
    }
    tensor2 operator * (const double &A)
    {
        tensor2 temp;
        temp.a11=this->a11*A;
        temp.a12=this->a12*A;
        temp.a13=this->a13*A;
        temp.a21=this->a21*A;
        temp.a22=this->a22*A;
        temp.a23=this->a23*A;
        temp.a31=this->a31*A;
        temp.a32=this->a32*A;
        temp.a33=this->a33*A;
        return temp;
    }    
        
    void operator = (const double A[6])
    {

        this->a11=A[0];
        this->a22=A[1];
        this->a33=A[2];
        this->a21=A[3];
        this->a12=A[3];
        this->a23=A[4];
        this->a32=A[4];
        this->a13=A[5];
        this->a31=A[5];

    }

};

 inline double trace_tensor2 (const tensor2 &A )
 {
   return (A.a11+A.a22+A.a33);
     
 };


 inline   tensor1  vectorproduct(const tensor1 &A,const tensor1 &B)
    {
        tensor1 temp;
        temp.a1=A.a2*B.a3-A.a3*B.a2;
        temp.a2=A.a3*B.a1-A.a1*B.a3;
        temp.a3=A.a1*B.a2-A.a2*B.a1;
        return temp;
    };

 inline   tensor2 kroneckerproduct(const tensor1 &A, const tensor1 &B)
    {
        tensor2 temp;
        temp.a11=A.a1*B.a1;
        temp.a12=A.a1*B.a2;
        temp.a13=A.a1*B.a3;
        temp.a21=A.a2*B.a1;
        temp.a22=A.a2*B.a2;
        temp.a23=A.a2*B.a3;
        temp.a31=A.a3*B.a1;
        temp.a32=A.a3*B.a2;
        temp.a33=A.a3*B.a3;
        return temp;
    };

inline    double doubleinnerproduct(const tensor2 &A, const tensor2 &B)
    {
      double temp;
      temp=A.a11*B.a11+A.a12*B.a21+A.a13*B.a31+A.a21*B.a12+A.a22*B.a22+A.a23*B.a32+A.a31*B.a13+A.a32*B.a23+A.a33*B.a33;
      return temp;
    };


 inline   double vonmises(const tensor2 &A)
    {
        double temp;
        temp= sqrt( ( (A.a11-A.a22)*(A.a11-A.a22) + (A.a22-A.a33)*(A.a22-A.a33) + (A.a33-A.a11)*(A.a33-A.a11) + 6.0*( A.a12*A.a12 + A.a23*A.a23+A.a31*A.a31) )*0.5) ;
        return temp;
    };



#endif
