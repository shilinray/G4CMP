#!/bin/bash

Al_z="100 200 300 400 500 600 700 800 900 1000"
Nb_z="20"
Al_bools="true false"
numsensors="1 2 3 4 5 6 7 8 9 10 11"

# Log file
logfile="loop_timing.log"

# Write overall start time
{
    echo "========================================"
    echo "loop started at: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "========================================"
} > "$logfile"

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

# Write overall end time
{
    echo "========================================"
    echo "loop ended at:   $(date '+%Y-%m-%d %H:%M:%S')"
    echo "========================================"
} >> "$logfile"
