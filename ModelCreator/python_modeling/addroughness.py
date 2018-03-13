#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 15:05:29 2017
!Add a semi-circle defect!
@author: xliang
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 17:56:43 2017

@author: xliang
"""

import os,re
limit=1
print("model contains %s defect(s)"%str(limit))

grainnumber=10000
width=2
height=1
chlength=0.005
radius=0.2550 # Attention 
str_rand = input("please enter the random seed(as the file name)\n")
#str_rand="5212"
str_name = "Model"+str_rand

#sys_com1 ="neper -T -n %s -dim 2 -morpho \"diameq:normal(14,1) inteval(12,16),sphericity:lognormal(0.145,0.03,1-x)\" -morphooptistop \"time=30,itermax=500000\" -regularization 1 -format geo,ori -domain \'square(%s,%s)\' -o %s -id %s" %(grainnumber,width,height,str_name,str_rand)

sys_com1="neper -T -n %s -dim 2 -format geo,ori -domain \'square(%s,%s)\' -o %s -id %s" %(grainnumber,width,height,str_name,str_rand)

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
filewrite=open("result_poly.geo",'w')
filewrite.write("""//+\n""")
filewrite.write("""SetFactory("OpenCASCADE");\n""")

for line in filepoly:
    if ("Physical" not in line):
        saveline(line,filewrite)

for i in range(0,plcount):
    i+=1
    filewrite.write("Plane Surface (%s) = {%s};\n"%(str(i),str(i)))

filewrite.write("For i In {1:%s} \n"%str(limit))
filewrite.write("   Disk(i + %s)={%s/(%s+1)*i,%s,0,%s};\n"%(str(plcount),str(width),str(limit),str(height),str(radius)))
filewrite.write("EndFor\n")
    
#for i in range(0,plcount-1):
    #i+=1
    #filewrite.write("BooleanDifference{Surface{%s};Delete;}{Surface{%s:%s};}\n"%(str(i),str(plcount+1),str(plcount+limit)))
    
filewrite.write("BooleanDifference{Surface{1:%s};Delete;}{Surface{%s:%s};Delete;}\n"%(str(plcount),str(plcount+1),str(plcount+limit)))
filewrite.write("Characteristic Length {1:newp-1} = %s;"%(chlength))

filepoly.close()
filewrite.close()

os.system("./ori2dat %s"%(str_name))
os.system("mv %s %s"%(str_name+".dat","grains_ori.dat"))
os.system("rm %s.ori"%(str_name))
