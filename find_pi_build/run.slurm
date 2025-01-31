#!/bin/bash
#SBATCH --job-name=find-pi-mpi
#SBATCH --partition=nu
#SBATCH --output=stdout-%x_%j
#SBATCH --error=stderr-%x_%j
#SBATCH --time=0-00:20:00
#SBATCH --nodes=1

module load kokkos/openmp
module load gcc/11.1.0
module load mpich/4.1.1
module load vtk/9.4.0
module load cmake/3.26.0

CXX=mpic++

$CXX -o find-pi-mpi findpi.cpp

mpiexec -np 4 ./find-pi-mpi 100000

rm find-pi-mi