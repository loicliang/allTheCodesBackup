 

import os,shutil

order="3"
filename="result_poly"

os.system("python3 addroughness.py")

cmdstring="gmsh -2 -order %s "%(order) + filename + ".geo"
os.system(cmdstring) 
os.system("./NewMSH2GEOF %s.msh"%(filename))
os.system("mv %s.geof %s"%(filename,"tt.geof"))
os.system("cp %s %s"%("tt.geof","model.geof"))



newfolder="0newcubmodel"
os.system("mkdir %s"%(newfolder))
os.system("cp %s %s"%("tt.geof",newfolder))
os.system("cp %s %s"%("model.geof",newfolder))
os.system("cp %s.msh %s/0.msh"%(filename,newfolder))
os.system("cp %s %s"%("repos/cub/CubElast.txt",newfolder+"/"))
os.system("cp %s %s"%("repos/cub/tt.inp",newfolder+"/"))
os.system("cp %s %s"%("repos/cub/newdodo.slurm",newfolder+"/"))
os.system("cp %s %s"%("repos/cub/fatOnly.slurm",newfolder+"/"))

os.system("cp %s %s"%("grains_ori.dat",newfolder))

newfolder+="/geoFiles"
os.system("mkdir %s"%(newfolder))
os.system("mv %s %s"%("*.geo",newfolder))


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
