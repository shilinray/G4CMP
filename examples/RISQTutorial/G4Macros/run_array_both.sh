#!/bin/bash

Al_z="100 200 300 400 500 600 700 800 900 1000"
Nb_z="20 40 60 80 100"
Al_bools="true false"

# First parameters: pAbsProbSideWallSi=0.01 pAbsProbPolishedWallSi=0.0
# Second parameters: pAbsProbSideWallSi=0.0 pAbsProbPolishedWallSi=0.0025

for alz in $Al_z
do
    for nbz in $Nb_z
    do
        for albool in $Al_bools
        do
            # Config 1
            echo $nbz $alz $albool 0.01 0.0
            sbatch ./thickness.sh $alz $nbz $albool 0.01 0.0
            
            # Config 2
            echo $nbz $alz $albool 0.0 0.0025
            sbatch ./thickness.sh $alz $nbz $albool 0.0 0.0025
        done
    done
done
