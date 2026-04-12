#!/bin/bash

Al_z="100 200 300 400 500 600 700 800 900 1000"
Nb_z="20 40 60 80 100"

for alz in $Al_z
do
    for nbz in $Nb_z
    do
        echo $nbz $alz
        sbatch ./thickness.sh $alz $nbz
    done
done

