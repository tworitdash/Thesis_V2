#!/bin/sh
#SBATCH --gres=gpu
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=END
#SBATCH -o fmincon_Opt_V2L_serv2.out

module use /opt/insy/modulefiles
module load cuda/9.0 cudnn/9.0-7.4.2.24 matlab/R2017b
srun matlab -r "fmincon_Opt_minxp_V2L_serv2; exit(0)"
