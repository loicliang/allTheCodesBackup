
import os,shutil

order="3"
filename="J2_defects_flexible-Final"

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
os.system("cp %s %s"%("dodoFEM_repos/isotropicElastic.txt",newfolder+"/"))
os.system("cp %s %s"%("dodoFEM_repos/fatigue.inp",newfolder+"/"))
os.system("cp %s %s"%("dodoFEM_repos/newdodo.slurm",newfolder+"/"))
os.system("cp %s %s"%("dodoFEM_repos/fatOnly.slurm",newfolder+"/"))
