#!/bin/bash

#SBATCH --partition=roma
#SBATCH --account=supercdms:default
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4096
#SBATCH --time=0-02:00:00
#SBATCH --gpus 0

# load SCDMS singularity image module

module load scdms/V05-02
module list

macro="pceStudy.mac"

echo "Macro is: '$macro'"

apptainer exec /sdf/group/supercdms/software/releases/cdmsfull_V05-02.sif \
  /sdf/home/s/shilin/mycode/G4CMP/examples/RISQTutorial/RISQTutorial-build/RISQTutorial \
  "$macro"
