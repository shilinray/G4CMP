#!/bin/bash

Al_z="100"
Nb_z="20"
Al_bools="false"
numsensors="1 2 5 10 20 30 40 50 60 70 80 90 100 125 150 175 200"


# First parameters:  pAbsProbSideWallSi=0.04  pAbsProbPolishedWallSi=0.0  (SW loss)
# Second parameters: pAbsProbSideWallSi=0.0   pAbsProbPolishedWallSi=0.01 (PF loss)

for alz in $Al_z
do
    for nbz in $Nb_z
    do
        for albool in $Al_bools
        do
            for ns in $numsensors
            do
                # Config 1
                echo "$nbz $alz $albool 0.04 0.0 $ns"
                sbatch ./numsensors.sh "$alz" "$nbz" "$albool" 0.04 0.0 "$ns"

                # Config 2
                echo "$nbz $alz $albool 0.0 0.01 $ns"
                sbatch ./numsensors.sh "$alz" "$nbz" "$albool" 0.0 0.01 "$ns"
            done
        done
    done
done

# macro chain, run_numsensors.sh then numsensors.sh then pceStudy.mac