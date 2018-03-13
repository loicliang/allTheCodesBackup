 
//+
SetFactory("OpenCASCADE");

// Set the parameters
//************************************************************
limit = DefineNumber[ 1 , Min 1, Max 36, Step 1,
  Name "Parameters/Defect Number" ];

X_lower_left_corner=0;  
Y_lower_left_corner=0; 
Z_lower_left_corner=0;
Width_Rectangle =   0.5;
Height_Rectangle =  0.5;

Radius_circle =  0.05;
Radius_pore   =  0.02;
//Main Program
//*********************************************************************************************************************************
Rectangle(1)={X_lower_left_corner,Y_lower_left_corner,Z_lower_left_corner,Width_Rectangle,Height_Rectangle}; //{x,y,z,width,height}

For i In {1:limit}
        
    Disk(i+1)={Width_Rectangle/(limit+1)*i,Height_Rectangle,0,Radius_circle};// {x,y,z,Radius,(for ellipse)}

EndFor

Disk(limit+2)={Width_Rectangle/2.0,Height_Rectangle/2.0,0 ,Radius_pore};// a circle in the center

BooleanDifference(limit+3) = { Surface{1};Delete; } { Surface{2:limit+1,limit+2 }; Delete; };

// Set the Characteristic Length
//**************************************************
Characteristic Length { 2,2*limit+4} = 0.1;
Characteristic Length { 1:2*limit+4 } = 0.002;

Characteristic Length { 2*limit+5 } = 0.005;
