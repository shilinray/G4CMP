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

if [ "$#" -lt 5 ]
then
    echo "Enter Al, Nb thicknesses, Al bool (true/false), pAbsProbSideWallSi, pAbsProbPolishedWallSi"
    # exit: quits the script immedietly
    exit
fi

# takes the first and second gives parameters and assigns them to variables
Al_z="$1"
Nb_z="$2"
Al_bool="$3"
pAbsProbSideWallSi="$4"
pAbsProbPolishedWallSi="$5"

if [ "$Al_bool" = "true" ]
then
    if [ "$pAbsProbPolishedWallSi" = "0.0" ]
    then
        run_dir="./260410_run/lQPD_Al_SW"
    elif [ "$pAbsProbSideWallSi" = "0.0" ]
    then
        run_dir="./260410_run/lQPD_Al_PF"
    else
        run_dir="./260410_run/lQPD_Al_wall_${pAbsProbSideWallSi}_polished_${pAbsProbPolishedWallSi}"
    fi
else
    if [ "$pAbsProbPolishedWallSi" = "0.0" ]
    then
        run_dir="./260410_run/lQPD_Al_Nb_SW"
    elif [ "$pAbsProbSideWallSi" = "0.0" ]
    then
        run_dir="./260410_run/lQPD_Al_Nb_PF"
    else
        run_dir="./260410_run/lQPD_Al_Nb_wall_${pAbsProbSideWallSi}_polished_${pAbsProbPolishedWallSi}"
    fi
fi 

tag="Al${Al_z}_Nb${Nb_z}"

echo "Al thickness: $Al_z"
echo "Nb thickness: $Nb_z"
echo "Al bool: $Al_bool"
echo "pAbsProbSideWallSi: $pAbsProbSideWallSi"
echo "pAbsProbPolishedWallSi: $pAbsProbPolishedWallSi"

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
sed -i "s/Al false/Al ${Al_bool}/g" $mName
sed -i "s/pAbsProbSideWallSi 0.0/pAbsProbSideWallSi ${pAbsProbSideWallSi}/g" $mName
sed -i "s/pAbsProbPolishedWallSi 0.0/pAbsProbPolishedWallSi ${pAbsProbPolishedWallSi}/g" $mName
sed -i "s|hitsFileName|hitsFileName ${fHits}|g" $mName
sed -i "s|primFileName|primFileName ${fPrimary}|g" $mName

apptainer exec /sdf/group/supercdms/software/releases/cdmsfull_V05-02.sif \
  /sdf/home/s/shilin/mycode/G4CMP/examples/RISQTutorial/RISQTutorial-build/RISQTutorial \
  "$mName"


