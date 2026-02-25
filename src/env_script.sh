#!/bin/bash
module load java/jdk-12_sequana
module load gcc/9.3_sequana
module load boost/1.87.0_gnu+openmpi-4.1.4_sequana
module load anaconda3/2024.02_sequana
module load masa-openmp/1.0.1.1024_sequana
#source activate /scratch/cenapadrjsd/rafael.terra2/conda-env
source /scratch/cenapadrjsd/COMPSs/3.3.3_gnu/compssenv
export PATH=$PATH:/scratch/cenapadrjsd/app/pa-star2/bin/
echo $PYTHONPATH