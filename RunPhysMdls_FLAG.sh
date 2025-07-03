#!/bin/env bash
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 144:0:0
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=10G
#SBATCH --account=kg98
#SBATCH --output=/fs04/kg98/stuarto/ComingUpShort/SLURM_OUTPUT/slurm-%j.out

MDL=$1
ADDMULT=0
LAW=0
TIMING=0

module load matlab/r2023b

CODEDIR="/fs04/kg98/stuarto/ComingUpShort"

matlab -nodisplay -nosplash -r "addpath(genpath(('${CODEDIR}'))); run_fitGNM_FLaG($MDL,$ADDMULT,$LAW); exit"

