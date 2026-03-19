#!/bin/bash

#SBATCH --partition=roma
#SBATCH --account=supercdms:default
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=21
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

run_dir="./260317_run/SQUAT"
tag="Al${Al_z}_Nb${Nb_z}"

echo "Al thickness: $Al_z"
echo "Nb thickness: $Nb_z"

mkdir -p "$run_dir"

fHits="$run_dir/Hits_${tag}.txt"
fPrimary="$run_dir/Primary_${tag}.txt"
mName="$run_dir/Run_${tag}.mac"

echo "Base Macro $BaseMacro"
echo "Macro File $mName"
echo "Hits File $fHits"
echo "Primary File $fPrimary"

# copies BaseMacro to a new name/file
cp $BaseMacro $mName

sed -i "s/filmThicknessAl 0/filmThicknessAl ${Al_z}/g" $mName
sed -i "s/filmThicknessNb 0/filmThicknessNb ${Nb_z}/g" $mName
sed -i "s|hitsFileName|hitsFileName ${fHits}|g" $mName
sed -i "s|primFileName|primFileName ${fPrimary}|g" $mName

apptainer exec /sdf/group/supercdms/software/releases/cdmsfull_V05-02.sif \
  /sdf/home/s/shilin/mycode/G4CMP/examples/RISQTutorial/RISQTutorial-build/RISQTutorial \
  "$mName"

mv ./260317_run/SQUAT /sdf/scratch/supercdms/sray/260317_run