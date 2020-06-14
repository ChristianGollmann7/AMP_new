#! /bin/bash

#SBATCH -p q_student
#SBATCH -N 1                 
#SBATCH -c 64   # use all 64 cores 
#SBATCH --cpu-freq=High
#SBATCH --time=5:00
#SBATCH --output=test_job.out


./Lock 9 60 100000 50

