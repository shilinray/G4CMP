#!/bin/bash

#SBATCH --partition=roma
#SBATCH --account=supercdms:default
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=4096
#SBATCH --time=0-02:00:00
#SBATCH --gpus 0

# load SCDMS singularity image module

module load scdms/V05-02
module list

BaseMacro="memsStudy.mac"
echo "Macro is: $BaseMacro"

if [ "$#" -lt 8 ]
then
    echo "Enter Al thickness, Nb thickness, Al bool (true/false), pAbsProbSideWallSi, pAbsProbPolishedWallSi, ix, iy, grid_size"
    exit
fi

Al_z="$1"
Nb_z="$2"
Al_bool="$3"
pAbsProbSideWallSi="$4"
pAbsProbPolishedWallSi="$5"
ix="$6"
iy="$7"
grid_size="$8"

# Compute x and y grid positions for grid_size x grid_size grid on 1cm x 1cm chip
# Spacing = 1/grid_size cm; first point at -0.5 + 0.5/grid_size cm
gps_x=$(echo "scale=4; -0.5 + ($ix + 0.5) / $grid_size" | bc -l)
gps_y=$(echo "scale=4; -0.5 + ($iy + 0.5) / $grid_size" | bc -l)

if [ "$Al_bool" = "true" ]
then
    if [ "$pAbsProbPolishedWallSi" = "0.0" ]
    then
        run_dir="./260519_run/MEMS_Al_SW"
    elif [ "$pAbsProbSideWallSi" = "0.0" ]
    then
        run_dir="./260519_run/MEMS_Al_PF"
    else
        run_dir="./260519_run/MEMS_Al_wall_${pAbsProbSideWallSi}_polished_${pAbsProbPolishedWallSi}"
    fi
else
    if [ "$pAbsProbPolishedWallSi" = "0.0" ]
    then
        run_dir="./260519_run/MEMS_Al_Nb_SW"
    elif [ "$pAbsProbSideWallSi" = "0.0" ]
    then
        run_dir="./260519_run/MEMS_Al_Nb_PF"
    else
        run_dir="./260519_run/MEMS_Al_Nb_wall_${pAbsProbSideWallSi}_polished_${pAbsProbPolishedWallSi}"
    fi
fi

tag="Al${Al_z}_Nb${Nb_z}_ix${ix}_iy${iy}"

echo "Al thickness: $Al_z"
echo "Nb thickness: $Nb_z"
echo "Al bool: $Al_bool"
echo "pAbsProbSideWallSi: $pAbsProbSideWallSi"
echo "pAbsProbPolishedWallSi: $pAbsProbPolishedWallSi"
echo "ix: $ix  iy: $iy"
echo "GPS position: ($gps_x, $gps_y) cm"

mkdir -p "$run_dir"

fHits="$run_dir/Hits_${tag}.txt"
fPrimary="$run_dir/Primary_${tag}.txt"
mName="$run_dir/Run_${tag}.mac"

echo "Base Macro $BaseMacro"
echo "Macro File $mName"
echo "Hits File $fHits"
echo "Primary File $fPrimary"

cp $BaseMacro $mName

sed -i "s/filmThicknessAl 0/filmThicknessAl ${Al_z}/g" $mName
sed -i "s/filmThicknessNb 0/filmThicknessNb ${Nb_z}/g" $mName
sed -i "s/Al false/Al ${Al_bool}/g" $mName
sed -i "s/pAbsProbSideWallSi 0.0/pAbsProbSideWallSi ${pAbsProbSideWallSi}/g" $mName
sed -i "s/pAbsProbPolishedWallSi 0.0/pAbsProbPolishedWallSi ${pAbsProbPolishedWallSi}/g" $mName
sed -i "s|GPS_X|${gps_x}|g" $mName
sed -i "s|GPS_Y|${gps_y}|g" $mName
sed -i "s|hitsFileName|hitsFileName ${fHits}|g" $mName
sed -i "s|primFileName|primFileName ${fPrimary}|g" $mName

apptainer exec /sdf/group/supercdms/software/releases/cdmsfull_V05-02.sif \
  /sdf/home/s/shilin/mycode/G4CMP/examples/RISQTutorial/RISQTutorial-build/RISQTutorial \
  "$mName"
