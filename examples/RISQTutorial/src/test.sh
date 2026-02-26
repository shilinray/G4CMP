#!/bin/bash

#SBATCH --partition=roma
#SBATCH --account=supercdms:default
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4096
#SBATCH --time=0-02:00:00
#SBATCH --gpus 0

echo "Hi there" > out.txt
