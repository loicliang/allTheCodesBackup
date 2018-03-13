 
//+
SetFactory("OpenCASCADE");

// Set the parameters
//************************************************************
limit = DefineNumber[ 1 , Min 1, Max 36, Step 1,
  Name "Parameters/Defect Number" ];

limit=1;  
X_lower_left_corner=0;  
Y_lower_left_corner=0; 
Z_lower_left_corner=0;
Width_Rectangle =   2;
Height_Rectangle =  1;

Radius_circle =  0.1825;

//Main Program
//*********************************************************************************************************************************
Rectangle(1)={X_lower_left_corner,Y_lower_left_corner,Z_lower_left_corner,Width_Rectangle,Height_Rectangle}; //{x,y,z,width,height}

For i In {1:limit}
        
    Disk(i+1)={Width_Rectangle/(limit+1)*i,Height_Rectangle,0,Radius_circle};// {x,y,z,Radius,(for ellipse)}

EndFor

BooleanDifference(limit+2) = { Surface{1};Delete; } { Surface{2:limit+1}; Delete; };

// Set the Characteristic Length
//**************************************************
Characteristic Length { 2,newp-1} = 0.2;
Characteristic Length { 1,newp-2} = 0.05;
Characteristic Length { 3:newp-3} = 0.001;
