#!/bin/bash

#SBATCH --time=0-00:10:00   # walltime limit (HH:MM:SS)
#SBATCH --partition=normal
#SBATCH --nodes=1   # number of nodes, KEEP AT 1 ON EULER, unless running NEB jobs
#SBATCH --ntasks-per-node=16   # 36 processor core(s) per node
#SBATCH --mem=175G   # maximum memory per node
#SBATCH --mail-user=vyn@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

ulimit -s unlimited

date > a
/home/dmpatel/zacros_4.0/build_mpi/zacros.x

date >> a

