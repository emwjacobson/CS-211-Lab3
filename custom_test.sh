#!/bin/bash
#SBATCH -o custom_test.o
#SBATCH -e custom_test.err
#SBATCH -J custom_test
#SBATCH -N 1
#SBATCH -t 00:05:00

module load mpich-3.2.1/gcc-4.8.5

mpirun -np 9 ./build/sieve3 10000000 &> output.txt
