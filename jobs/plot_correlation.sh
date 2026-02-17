#!/bin/bash
#SBATCH -J plot_corr  # run's name
#SBATCH -N 1                    # request 1 node
#SBATCH -c 1                  # request 1 cpu per task
#SBATCH --mem=20GB                 # request 10GB
#SBATCH -t 01:00:00                # request 12 hours walltime
#SBATCH -o Out_plots.txt                 # output file name
#SBATCH -e Err_plots.txt                 # error file name
#SBATCH --mail-type=BEGIN,END,FAIL # send me a mail at beginning and end of the job
#SBATCH --mail-user=laurent.magri-stella@lam.fr


julia /net/GECO/nas12c/euclid/TestWorkFolder/correlation_plots.jl

