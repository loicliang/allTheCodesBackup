
import os,shutil

order="3"
filename="J2_defects_flexible-Final"

cmdstring="gmsh -2 -order %s "%(order) + filename + ".geo"
os.system(cmdstring) 
os.system("./NewMSH2GEOF %s.msh"%(filename))
os.system("mv %s.geof %s"%(filename,"tt.geof"))
os.system("cp %s %s"%("tt.geof","model.geof"))


newfolder="0newisomodel"
os.system("mkdir %s"%(newfolder))
os.system("mv %s %s"%("tt.geof",newfolder))
os.system("mv %s %s"%("model.geof",newfolder))
os.system("mv %s.msh %s/0iso.msh"%(filename,newfolder))
os.system("cp %s %s"%("dodoFEM_repos/iso/isotropicElastic.txt",newfolder+"/"))
os.system("cp %s %s"%("dodoFEM_repos/iso/fatigue.inp",newfolder+"/"))
os.system("cp %s %s"%("dodoFEM_repos/iso/newdodo.slurm",newfolder+"/"))
os.system("cp %s %s"%("dodoFEM_repos/iso/fatOnly.slurm",newfolder+"/"))
