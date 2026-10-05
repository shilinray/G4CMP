#!/bin/bash

Al_z="100"
Nb_z="20"
Al_bools="false"
numsensors="1"
# 100ueV, 200ueV, 300ueV, 600ueV, 1meV, 2meV, 3meV, 4meV, 5meV, 10meV (in eV)
energies="0.0001 0.0002 0.0003 0.0006 0.001 0.002 0.003 0.004 0.005 0.01"

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
                    echo "$nbz $alz $albool 0.01 0.0 $ns $e"
                    sbatch ./energy.sh "$alz" "$nbz" "$albool" 0.01 0.0 "$ns" "$e"

                    # Config 2 (PF loss)
                    echo "$nbz $alz $albool 0.0 0.0025 $ns $e"
                    sbatch ./energy.sh "$alz" "$nbz" "$albool" 0.0 0.0025 "$ns" "$e"
                done
            done
        done
    done
done

# macro chain, run_energy.sh then energy.sh then pceStudy.mac
