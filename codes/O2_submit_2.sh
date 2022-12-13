#!/bin/bash
#SBATCH -c 4                               # Request one core
#SBATCH --job-name=HST508_final_project
#SBATCH -t 0-6:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=16GB                         # Memory total in GB (for all cores)
#SBATCH -o HST508_final_project%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e HST508_final_project%j.err                 # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-user=marie_zhang@hms.harvard.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load miniconda3/4.10.3

source activate test

python3 hic.py --run_analysis1 ${1} --run_analysis2 ${2} --run_analysis3 ${3} --run_analysis4 ${4}

conda deactivate

# sbatch O2_submit_2.sh False False False True