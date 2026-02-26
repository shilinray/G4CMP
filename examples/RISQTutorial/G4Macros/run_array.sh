#!/bin/bash

dead_z="20 40 60 80 100"
live_z="100 200 300 400 500 600"

for dz in $dead_z
do
    for lz in $live_z
    do
        echo $dz $lz
        sbatch ./thickness.sh $dz $lz
    done
done


