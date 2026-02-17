#!/bin/bash
#SBATCH -J comp_flexion_corr_50arcm  # run's name
#SBATCH -N 1                    # request 1 node
#SBATCH -c 20                  # request 1 cpu per task
#SBATCH --mem=240GB                 # request 10GB
#SBATCH -t 08:00:00                # request 12 hours walltime
#SBATCH -o Out.txt                 # output file name
#SBATCH -e Err.txt                 # error file name
#SBATCH --mail-type=BEGIN,END,FAIL # send me a mail at beginning and end of the job
#SBATCH --mail-user=laurent.magri-stella@lam.fr

export JULIA_NUM_THREADS=20

julia /net/GECO/nas12c/euclid/TestWorkFolder/correlation_function.jl

