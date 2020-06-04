#! /bin/bash

#SBATCH -p q_student
#SBATCH -N 1                 
#SBATCH -c 32   # use all 32 cores 
#SBATCH --cpu-freq=High
#SBATCH --time=5:00
#SBATCH --output=test_job.out


./Lock 9 2 100 20

