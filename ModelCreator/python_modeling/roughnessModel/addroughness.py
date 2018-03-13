#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 15:05:29 2017

@author: xliang
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 17:56:43 2017

@author: xliang
"""

import os,re

grainnumber=6000
width=2.0
height=0.6
chlength=0.005
str_rand = input("please enter the random seed(as the file name)\n")
str_name = "Model"+str_rand
sys_com1 ="neper -T -n %s -dim 2 -morpho \"diameq:normal(15,1) inteval(12,18),sphericity:lognormal(0.145,0.03,1-x)\" -morphooptistop \"time=3000,itermax=500000\" -regularization 1 -format geo,ori -domain \'square(%s,%s)\' -o %s -id %s" %(grainnumber,width,height,str_name,str_rand)
os.system(sys_com1)


pointline=re.compile(r'Point\s*\((.*)\)\s*\=\s*\{(.*)\,(.*)\,(.*)') 
lineline=re.compile(r'Line\s*\((.*)\)\s*\=\s*\{(.*)\,(.*)\}')
surfaceline=re.compile(r'Line Loop\s*\((.*)\)\s*\=\s*\{(.*)\}')
ptcount=0
lncount=0
plcount=0


def saveline(a_string,fout):
    if(pointline.match(a_string)):
        global ptcount
        ptcount+=1
        fout.write(a_string)
    elif(lineline.match(a_string)):
        global lncount
        lncount+=1
        fout.write(a_string)
    elif(surfaceline.match(a_string)):
        global plcount
        plcount+=1
        fout.write(a_string)
    else:
        pass

pointdict={}
linedict={}
planedict={}

def addrough(a_string,fout):
    if(pointline.match(a_string)):
        global ptcount
        ptcount+=1
        pointdict[(pointline.match(a_string).group(1))]=ptcount
        a_string=re.sub(r'\((.*)\)','('+str(ptcount)+')',a_string)
        fout.write(a_string)
    elif(lineline.match(a_string)):
        global lncount
        lncount+=1
        linedict[(lineline.match(a_string).group(1))]=lncount
        seq=re.search(r'\{(.*)\}',a_string)
        seq=re.split(r'[\s\,]+',seq.group(1))
        outline="Line(%s) = {%s,%s};\n"%(str(lncount),str(pointdict[seq[0]]),str(pointdict[seq[1]]))
        fout.write(outline)
    elif(surfaceline.match(a_string)):
        global plcount
        plcount+=1
        planedict[(surfaceline.match(a_string).group(1))]=plcount
        seq=re.search(r'\{(.*)\}',a_string)
        seq=re.split(r'[\s\,]+',seq.group(1))
        outline="Line Loop(%s) = {"%(plcount)
        for i in seq:
            outline+=str(linedict[i])
            outline+=','
        outline=outline[:-1]
        outline+="};\n"
        fout.write(outline)
    else:
        pass
    
    
    
                
filepoly=open(str_name+".geo")
filerough=open("rough.geo")
filewrite=open("result_poly.geo",'w')
filewrite.write("""//+\n""")
filewrite.write("""SetFactory("OpenCASCADE");\n""")

for line in filepoly:
    if ("Physical" not in line):
        saveline(line,filewrite)
plpoly=plcount
for line in filerough:
    if ("Physical" not in line):
        addrough(line,filewrite)

print(ptcount,lncount,plcount)

for i in range(0,plcount):
    i+=1
    filewrite.write("Plane Surface (%s) = {%s};\n"%(str(i),str(i)))

filewrite.write("Disk(%s)={%s,%s,0,0.05};\n"%(str(plcount+1),str(width/2.0),str(height/2.0)))
plcount+=1
        
filewrite.write("BooleanDifference{Surface{1:%s};Delete;}{Surface{%s:%s};Delete;}\n"%(str(plpoly),str(plpoly+1),str(plcount)))
filewrite.write("Characteristic Length {1:newp-1} = %s;"%(chlength))



filerough.close()
filepoly.close()
filewrite.close()

os.system("./ori2dat %s"%(str_name))
os.system("mv %s %s"%(str_name+".dat","grains_ori.dat"))
os.system("rm %s.ori"%(str_name))
