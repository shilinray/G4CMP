#!/bin/bash

#SBATCH --partition=roma
#SBATCH --account=supercdms:default
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4096
#SBATCH --time=0-02:00:00
#SBATCH --gpus 0

# load SCDMS singularity image module

module load scdms/V05-02
module list

BaseMacro="pceStudy.mac"
echo "Macro is: '$macro'"

if [ "$#" -lt 2 ]
then
    echo "Enter dead and live thicknesses"
    # exit: quits the script immedietly
    exit
fi

# takes the first and second gives parameters and assigns them to variables
dead_z = $1
live_z = $2

echo "Dead Metal thickness: "$dead_z
echo "Live Metal thickness: "$live_z

# creates output files names and macro file name
fHits="Hits_dz"${dead_z}"_lz"${live_z}".txt"
fPrimary="Primary_dz"${dead_z}"_lz"${live_z}".txt"
mName="Run_dz"${dead_z}"_lz"${live_z}".txt"

# prints out a summary
echo "Base Macro "$BaseMacro
echo "Macro File"$mName
echo "Hits File"$fHits
echo "Primary File"$fPrimary

# copies BaseMacro to a new name/file
cp $BaseMacro $mName

sed -i

apptainer exec /sdf/group/supercdms/software/releases/cdmsfull_V05-02.sif \
  /sdf/home/s/shilin/mycode/G4CMP/examples/RISQTutorial/RISQTutorial-build/RISQTutorial \
  "$macro"
