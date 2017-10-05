 
//+
SetFactory("OpenCASCADE");

// Set the parameters
//************************************************************
limit = DefineNumber[ 1 , Min 1, Max 36, Step 1,
  Name "Parameters/Defect Number" ];

X_lower_left_corner=0;  
Y_lower_left_corner=0; 
Z_lower_left_corner=0;
Width_Rectangle =   4;
Height_Rectangle =  4;

Radius_circle =  0.05;

//Main Program
//*********************************************************************************************************************************
Rectangle(1)={X_lower_left_corner,Y_lower_left_corner,Z_lower_left_corner,Width_Rectangle,Height_Rectangle}; //{x,y,z,width,height}

For i In {1:limit}
        
    Disk(i+1)={Width_Rectangle/(limit+1)*i,Height_Rectangle,0,Radius_circle};// {x,y,z,Radius,(for ellipse)}

EndFor

BooleanDifference(limit+2) = { Surface{1};Delete; } { Surface{2:limit+1}; Delete; };

// Set the Characteristic Length
//**************************************************
Characteristic Length { 2,2*limit+4} = 1;
Characteristic Length { 1,2*limit+3} = 0.2;
Characteristic Length { 1:2*limit+4 } = 0.001;
