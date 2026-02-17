#!/bin/bash
#SBATCH -J rebuild_lensing_maps  # run's name
#SBATCH -N 1                    # request node
#SBATCH -c 1                # request cpu per task
#SBATCH --mem=240GB                 # request 10GB
#SBATCH -t 1-00:00:00                # request 12 hours walltime
#SBATCH -o Out_maps.txt                 # output file name
#SBATCH -e Err_maps.txt                 # error file name
#SBATCH --mail-type=BEGIN,END,FAIL # send me a mail at beginning and end of the job
#SBATCH --mail-user=laurent.magri-stella@lam.fr

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK 

julia /net/GECO/nas12c/euclid/TestWorkFolder/rebuild_maps.jl

