//+
SetFactory("OpenCASCADE");
 
X_lower_left_corner=0;  
Y_lower_left_corner=0; 
Z_lower_left_corner=0;
Width_Rectangle =   4;
Height_Rectangle =  4;


Rectangle(1)={X_lower_left_corner,Y_lower_left_corner,Z_lower_left_corner,Width_Rectangle,Height_Rectangle}; //{x,y,z,width,height}

// Set the Characteristic Length
//**************************************************
Characteristic Length { 1,2} = 0.2;
Characteristic Length { 1:4 } = 0.01;
