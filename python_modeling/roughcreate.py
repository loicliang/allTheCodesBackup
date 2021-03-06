#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 17:25:03 2017

@author: xliang
"""
import math
pi=3.1415926535897932

file=open("rough.geo",'w')


pts=100

for i in range(0,pts+1):
    x=10.0/pts*i
    y=math.sin(i*pi/10)*0.05+0.5
    file.write("Point(%s) = {%s,%s,0.0};\n"%(i+1,x,y))
file.write("Point(%s) = {10,10,0};\n"%(pts+2))    
file.write("Point(%s) = {0,10,0};\n"%(pts+3))  
for i in range(0,pts):
    file.write("Line(%s) = {%s,%s};\n"%(i+1,i+1,i+2))
file.write("Line(%s) = {%s,%s};\n"%(pts+1,pts+1,pts+2))
file.write("Line(%s) = {%s,%s};\n"%(pts+2,pts+2,pts+3))
file.write("Line(%s) = {%s,%s};\n"%(pts+3,pts+3,1))

outline="Line Loop (1) ={"
for i in range(1,pts+4):
    outline+=str(i)+','
outline=outline[:-1]
outline+="};"
file.write(outline)
file.close()