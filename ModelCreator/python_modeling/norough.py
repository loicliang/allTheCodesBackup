 

import os,re

grainnumber=10000
width=2
height=1
chlength=0.01
str_rand = input("please enter the random seed(as the file name)\n")
#str_rand="5212"
str_name = "Model"+str_rand
sys_com1 ="neper -T -n %s -dim 2 -morpho \"diameq:normal(14,1) inteval(12,16),sphericity:lognormal(0.145,0.03,1-x)\" -morphooptistop \"time=30,itermax=500000\" -regularization 1 -format geo,ori -domain \'square(%s,%s)\' -o %s -id %s" %(grainnumber,width,height,str_name,str_rand)
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

filewrite.write("Characteristic Length {1:newp-1} = %s;"%(chlength))

filepoly.close()
filewrite.close()

os.system("./ori2dat %s"%(str_name))
os.system("mv %s %s"%(str_name+".dat","grains_ori.dat"))
os.system("rm %s.ori"%(str_name))

order="3"
filename="result_poly"


cmdstring="gmsh -2 -order %s "%(order) + filename + ".geo"
os.system(cmdstring) 
os.system("./NewMSH2GEOF %s.msh"%(filename))
os.system("mv %s.geof %s"%(filename,"tt.geof"))
os.system("cp %s %s"%("tt.geof","model.geof"))



newfolder="0newaggmodel"
os.system("mkdir %s"%(newfolder))
os.system("mv %s %s"%("tt.geof",newfolder))
os.system("mv %s %s"%("model.geof",newfolder))
os.system("mv %s.msh %s/0.msh"%(filename,newfolder))
os.system("cp %s %s"%("repos/cp/316L_CP-mod.txt",newfolder+"/"))
os.system("cp %s %s"%("repos/cp/tt.inp",newfolder+"/"))
os.system("cp %s %s"%("repos/cp/newdodo.slurm",newfolder+"/"))
os.system("cp %s %s"%("repos/cp/fatOnly.slurm",newfolder+"/"))

os.system("mv %s %s"%("grains_ori.dat",newfolder))

newfolder+="/geoFiles"
os.system("mkdir %s"%(newfolder))
os.system("mv %s %s"%("*.geo",newfolder))
