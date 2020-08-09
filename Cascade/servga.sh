#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=END
#SBATCH -o Ga_Opt_RL_V2L.out

module use /opt/insy/modulefiles
module load cuda/9.0 cudnn/9.0-7.4.2.24 matlab/R2017b
srun matlab -r "Ga_Opt_RL_V2L; exit(0)"
