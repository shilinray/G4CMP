#!/bin/bash

Al_z="100 300 500 700 900"
Nb_z="20"
Al_bools="true false"
numsensors="1 5 10 20 40 80 100"


# First parameters: pAbsProbSideWallSi=0.01 pAbsProbPolishedWallSi=0.0
# Second parameters: pAbsProbSideWallSi=0.0 pAbsProbPolishedWallSi=0.0025

for alz in $Al_z
do
    for nbz in $Nb_z
    do
        for albool in $Al_bools
        do
            for ns in $numsensors
            do
                # Config 1
                echo "$nbz $alz $albool 0.01 0.0 $ns"
                sbatch ./numsensors.sh "$alz" "$nbz" "$albool" 0.01 0.0 "$ns"

                # Config 2
                echo "$nbz $alz $albool 0.0 0.0025 $ns"
                sbatch ./numsensors.sh "$alz" "$nbz" "$albool" 0.0 0.0025 "$ns"
            done
        done
    done
done