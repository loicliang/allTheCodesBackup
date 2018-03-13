#! /bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks-per-node 6
#SBATCH -p work
#SBATCH --mem=40G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=2016-2060@ensam.eu 
source /opt/share/intel/bin/ifortvars.sh intel64
module load mpi/mpich-intel
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
rm -f tt.out
srun --mpi=pmi2 /home/xliang/INSTALL6/bin/dodoFEM tt
/home/xliang/postTreatment/trunk/fatigue tt 4 5
/home/xliang/INSTALL6/bin/dodo2vtk tt

