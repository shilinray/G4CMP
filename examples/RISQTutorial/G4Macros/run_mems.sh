#!/bin/bash

Al_z="30"
Nb_z="30"
Al_bools="false"
grid_size=20  # number of points along each axis (grid_size x grid_size total)

grid_indices=$(seq 0 $((grid_size - 1)))

# First parameters:  pAbsProbSideWallSi=0.04  pAbsProbPolishedWallSi=0.0  (SW loss)
# Second parameters: pAbsProbSideWallSi=0.0   pAbsProbPolishedWallSi=0.01 (PF loss)

for albool in $Al_bools
do
    for ix in $grid_indices
    do
        for iy in $grid_indices
        do
            # Config 1: SW loss
            echo "$Al_z $Nb_z $albool 0.04 0.0 $ix $iy $grid_size"
            sbatch ./mems.sh "$Al_z" "$Nb_z" "$albool" 0.01 0.0 "$ix" "$iy" "$grid_size"

            # Config 2: PF loss
            echo "$Al_z $Nb_z $albool 0.0 0.01 $ix $iy $grid_size"
            sbatch ./mems.sh "$Al_z" "$Nb_z" "$albool" 0.0 0.0025 "$ix" "$iy" "$grid_size"
        done
    done
done


# macro chain, run_mems.sh then mems.sh then memsStudy.mac