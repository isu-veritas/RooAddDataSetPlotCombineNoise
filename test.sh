#!/bin/bash
#Submit this script with: sbatch thefilename
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-user=achrmy@iastate.edu
#SBATCH --mail-type=END

../bin/AddRooDataSet /work/LAS/amandajw-lab/users/achrmy/macros/20deg/itm/grisudetATM21ITMzen20.list 34 > log/grisudetATM21ITMzen20.log 2>&1
