#!/bin/bash

# Macro chain for the circular chip: run_mems_circle.sh then mems_circle.sh then memsStudy.mac
# Phonons are deposited on a grid_size x grid_size grid spanning the square that circumscribes
# the chip (side = 2 x chip_radius).

Al_z="30"
Nb_z="70"
Al_bools="false"
grid_size=50  # number of points along each axis (grid_size x grid_size total)
chip_radius=3.81  # cm (1.5 inch)
skip_outside_circle=false  # grid points outside the circular chip are not in Si, so don't submit them

grid_indices=$(seq 0 $((grid_size - 1)))

# SW loss:  pAbsProbSideWallSi=0.0088  pAbsProbPolishedWallSi=0.0
# PF loss:  pAbsProbSideWallSi=0.0   pAbsProbPolishedWallSi=0.000265

for albool in $Al_bools
do
    for ix in $grid_indices
    do
        for iy in $grid_indices
        do
            if [ "$skip_outside_circle" = "true" ]
            then
                inside=$(echo "scale=6; x = -$chip_radius + ($ix + 0.5) * 2 * $chip_radius / $grid_size; y = -$chip_radius + ($iy + 0.5) * 2 * $chip_radius / $grid_size; if (x*x + y*y <= $chip_radius*$chip_radius) 1 else 0" | bc -l)
                if [ "$inside" != "1" ]
                then
                    continue
                fi
            fi

            # Config 1: SW loss
            echo "$Al_z $Nb_z $albool 0.0088 0.0 $ix $iy $grid_size $chip_radius"
            sbatch ./mems_circle.sh "$Al_z" "$Nb_z" "$albool" 0.0088 0.0 "$ix" "$iy" "$grid_size" "$chip_radius"

            # Config 2: PF loss
            echo "$Al_z $Nb_z $albool 0.0 0.000265 $ix $iy $grid_size $chip_radius"
            sbatch ./mems_circle.sh "$Al_z" "$Nb_z" "$albool" 0.0 0.000265 "$ix" "$iy" "$grid_size" "$chip_radius"

            # Config 3: Same loss (side wall and polished face)
            echo "$Al_z $Nb_z $albool 0.000265 0.000265 $ix $iy $grid_size $chip_radius"
            sbatch ./mems_circle.sh "$Al_z" "$Nb_z" "$albool" 0.000265 0.000265 "$ix" "$iy" "$grid_size" "$chip_radius"
        done
    done
done
