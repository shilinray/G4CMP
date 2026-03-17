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
echo "Macro is: $BaseMacro"

if [ "$#" -lt 2 ]
then
    echo "Enter Al and Nb thicknesses"
    # exit: quits the script immedietly
    exit
fi

# takes the first and second gives parameters and assigns them to variables
Al_z="$1"
Nb_z="$2"

echo "Al thickness: $Al_z"
echo "Nb thickness: $Nb_z"

# creates output files names and macro file name
fHits="/sdf/scratch/supercdms/sray/260317_run/Hits_Al${Al_z}_Nb${Nb_z}.txt"
fPrimary="/sdf/scratch/supercdms/sray/260317_run/Primary_Al${Al_z}_Nb${Nb_z}.txt"
mName="Run_Al${Al_z}_Nb${Nb_z}.mac"

# prints out a summary
echo "Base Macro $BaseMacro"
echo "Macro File $mName"
echo "Hits File $fHits"
echo "Primary File $fPrimary"

# copies BaseMacro to a new name/file
cp $BaseMacro $mName

sed -i "s/filmThicknessAl 0/filmThicknessAl ${Al_z}/g" $mName
sed -i "s/filmThicknessNb 0/filmThicknessNb ${Nb_z}/g" $mName
sed -i "s|hitsFileName|hitsFileName '${fHits}'|g" $mName
sed -i "s|primFileName|primFileName '${fPrimary}'|g" $mName

apptainer exec /sdf/group/supercdms/software/releases/cdmsfull_V05-02.sif \
  /sdf/home/s/shilin/mycode/G4CMP/examples/RISQTutorial/RISQTutorial-build/RISQTutorial \
  "$mName"
