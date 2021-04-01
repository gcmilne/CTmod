#!/bin/bash

#SBATCH --output=/storage/users/gmilne/test/error/%x.out
#SBATCH --error=/storage/users/gmilne/test/error/%x.err
#SBATCH --time=005:00:00
#SBATCH --mem=5G

module load apps/R-3.6.3.tcl

cd /storage/users/gmilne/test 
R CMD BATCH --no-save /storage/users/gmilne/test/fitting.R
