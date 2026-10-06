#!/bin/bash

Al_z="40"
Nb_z="40"
Al_bools="false"
numsensors="1"
# (in eV)
energies="0.0002 0.0003 0.0004 0.0005 0.0006 0.0008 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.0075 0.01"

#1x1 chip
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
                for e in $energies
                do
                    # Config 1 (SW loss)
                    echo "$nbz $alz $albool 0.04 0.0 $ns $e"
                    sbatch ./energy.sh "$alz" "$nbz" "$albool" 0.04 0.0 "$ns" "$e"

                    # Config 2 (PF loss)
                    echo "$nbz $alz $albool 0.0 0.01 $ns $e"
                    sbatch ./energy.sh "$alz" "$nbz" "$albool" 0.0 0.01 "$ns" "$e"
                done
            done
        done
    done
done

# macro chain, run_energy.sh then energy.sh then pceStudy.mac
