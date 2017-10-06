/*
 * file surface_function.h
 * brief convert mesh to graph
 * author lxy
 * date 20170830 
*/

#ifndef _SURFACE_FUNCTION_H
#define _SURFACE_FUNCTION_H

#include <math.h>

#define err_lxy 0.0001 //compare the coordinates


#define X_center 2
#define Y_center 4
#define Z_center 0

/*
#define X_center 0.25
#define Y_center 0.45
#define Z_center 0
*/
#define DEPTH_SEARCH 5   // find the neighbour of the center
#define ZONE_RADIUS 0.02 // size of circle


#define MAXIUM_CYCLE 1   // the same as our loading

double distance_fun(double x, double y, double z)
{
    double res;
    res= sqrt ( (x-X_center)*(x-X_center)+(y-Y_center)*(y-Y_center)+(z-Z_center)*(z-Z_center) );
    return res;
    
};

double surface_function(double x) //1 defect, middle, size 50um
{
    double y;
/*   
    if( (x<=1.95) or (x>=2.05) )
    {
        y=4.0;
    }
    else
    {
        //y=4;
        y=4.0 - sqrt( ( 0.05*0.05 )-(x-2.0)*(x-2.0) );
    }
 */  
 /*
    if( (x<=0.2) or (x>=0.3) )
    {
        y=0.5;
    }
    else
    {
        y=0.5 - sqrt( ( 0.05*0.05 )-(x-0.25)*(x-0.25) );
    }
    */
 
 y=4.0;
    return y;
};



#endif
