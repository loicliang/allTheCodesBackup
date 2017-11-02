#not good No Fset in MSH2GEOF because not in OpenCascade.....
import os,shutil

order="3"
filename="result_poly"

grainnumber=6000
width=4.0
height=0.6
chlength=0.001
str_rand = input("please enter the random seed(as the file name)\n")
str_name = "Model"+str_rand
sys_com1 ="neper -T -n %s -dim 2 -morpho \"diameq:normal(15,1) inteval(12,18),sphericity:lognormal(0.145,0.03,1-x)\" -morphooptistop \"time=3000,itermax=500000\" -regularization 1 -format geo,ori -domain \'square(%s,%s)\' -o %s -id %s" %(grainnumber,width,height,str_name,str_rand)
os.system(sys_com1)

os.system("./ori2dat %s"%(str_name))
os.system("mv %s %s"%(str_name+".dat","grains_ori.dat"))
os.system("rm %s.ori"%(str_name))
os.system("mv %s.geo %s.geo"%(str_name,filename))

cmdstring="gmsh -2 -order %s "%(order) + filename + ".geo"
os.system(cmdstring) 
os.system("./NewMSH2GEOF %s.msh"%(filename))
os.system("mv %s.geof %s"%(filename,"tt.geof"))
os.system("cp %s %s"%("tt.geof","model.geof"))



newfolder="0newmodel"
os.system("mkdir %s"%(newfolder))
os.system("mv %s %s"%("tt.geof",newfolder))
os.system("mv %s %s"%("model.geof",newfolder))
os.system("mv %s.msh %s/0.msh"%(filename,newfolder))
os.system("cp %s %s"%("repos/316L_CP-mod.txt",newfolder+"/"))
os.system("cp %s %s"%("repos/tt.inp",newfolder+"/"))
os.system("cp %s %s"%("repos/newdodo.slurm",newfolder+"/"))
os.system("cp %s %s"%("repos/fatOnly.slurm",newfolder+"/"))

os.system("mv %s %s"%("grains_ori.dat",newfolder))

newfolder+="/geoFiles"
os.system("mkdir %s"%(newfolder))
os.system("mv %s %s"%("*.geo",newfolder))
