#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --job-name=es_all      	    # Job name
#SBATCH -N 3
#SBATCH --mem-per-cpu=20G           # Memory needed per CPU
#SBATCH --output=wes.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=wes.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE            # Type of email notification (BEGIN, END, FAIL, ALL)

# Automatically create a processing description
# bcbio_nextgen.py -w template ../config/wes-template.yaml name.csv WES/

bcbio_nextgen.py ../config/name.yaml -n 120 -t ipython -s slurm -q medium --tag "wes"

